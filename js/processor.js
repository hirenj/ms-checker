var fs      = require('fs'),
    sqlite3 = require('sqlite3');

var nconf = require('nconf');


Math.median = require('../js/math-median').median;

var uniprot_meta = require('../js/uniprot.js');

var promisify_sqlite = require('../js/sqlite-promise');

var quantitative = require('../js/quantitative');
var ambiguous = require('../js/ambiguous');
var hexnac_hcd = require('../js/hexnac_hcd');
var fragmentation = require('../js/fragmentions');
var metadata = require('../js/metadata');

var peptide = require('../js/peptide');

var spectra = require('../js/spectrum');

ambiguous.conf = nconf;
hexnac_hcd.conf = nconf;

[quantitative,ambiguous,hexnac_hcd,fragmentation,peptide].forEach(function(module) {
    module.on('task',function(desc) {
        console.log(module.constructor.name,desc);
        this.addListener('progress',function(percentage) {
            console.log(module.constructor.name,desc, (100*percentage).toFixed(0));
            if (percentage == 1) {
                module.removeListener('progress',arguments.callee);
            }
        });
    });
});


var partial = function(fn) {
    var args = Array.prototype.slice.call(arguments, 1);
    return function() { return fn.apply(this, args.concat(Array.prototype.slice.call(arguments, 0))); };
};

var onlyUnique = function(value, index, self) {
    return self.indexOf(value) === index;
};

var open_db = function(filename) {
    return new Promise(function(resolve,reject) {
        var db = new sqlite3.Database(filename,sqlite3.OPEN_READONLY,function(err) {
            if (err) {
                reject(err);
                return;
            }
            resolve(db);
        });
    });
};


var process_data = function(filename) {
    var global_datablock = {'data' : {}, 'metadata' : {}};
    return open_db(filename).then(function(db) {
        promisify_sqlite(db);
        metadata.get_metadata(db).then(function(meta) {
            global_datablock.metadata = meta;
        });
        var quant_promise = quantitative.init_caches(db).then(function() {
             return quantitative.retrieve_quantified_peptides(db).then(partial(quantitative.check_quantified_peptides,db));
        });
        var ambig_promise = ambiguous.retrieve_peptides(db).then(function(peps) { return peps? peps : []; });

        var processing_promise = Promise.all([quant_promise,ambig_promise]).then(function(all_peps) {
            var merged = Array.prototype.concat.apply([], all_peps);
            return merged;
        }).then(peptide.filter_ambiguous_spectra)
          .then(peptide.produce_peptide_scores_and_cleanup.bind(peptide,db))
          .then(hexnac_hcd.guess_hexnac.bind(hexnac_hcd,db))
          .then(partial(fragmentation.validate_peptide_coverage,db));

        return processing_promise.then(function(peps) {
            return uniprot_meta.init().then(function() {
                global_datablock.data = peptide.combine(peps);
                quantitative.clear_caches();
                return global_datablock;
            });
        });
    });
};

var process_files = function(files) {
    var result = Promise.resolve(true);
    var blocks = [];
    return result.then(function() {
        var self_func = arguments.callee;
        var file = files.shift();
        if (! file) {
            return blocks;
        }
        return process_data(file).then(function(block) {
            blocks.push(block);
            return self_func();
        });
    });
};

var combine = function(blocks,sources) {
    var result = { data : {} , metadata : [] };
    blocks.forEach(function(block) {
        var source = (sources || []).shift();

        Object.keys(block.data).forEach(function(prot) {
            if ( ! result.data[prot] ) {
                result.data[prot] = [];
            }
            if (source) {
                block.data[prot].forEach(function(pep) {
                    pep.source = source;
                });
            }
            result.data[prot] = result.data[prot].concat(block.data[prot]);
        });
        result.metadata.push(block.metadata);
    });
    return result;
};

module.exports = exports = { process: process_files, combine: combine };

