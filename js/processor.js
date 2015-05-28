var fs      = require('fs'),
    sqlite3 = require('sqlite3');

var nconf = require('nconf');
var ProgressBar = require('progress');

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

var ProgressBar = require('progress');

ambiguous.conf = nconf;
hexnac_hcd.conf = nconf;

var progress_bars = {};

function Multibar(stream) {
  this.stream     = stream || process.stderr;
  this.cursor     = 0;
  this.bars       = [];
  this.terminates = 0;
};

Multibar.prototype = {
  newBar: function(schema, options) {
    options.stream = this.stream;
    var bar = new ProgressBar(schema, options);
    this.bars.push(bar);
    var index = this.bars.length - 1;

    // alloc line
    this.move(index);
    this.stream.write('\n');
    this.cursor ++;

    // replace original
    var self  = this;
    bar.otick = bar.tick;
    bar.oterminate = bar.terminate;
    bar.tick = function(value, options) {
      self.tick(index, value, options);
    }
    bar.terminate = function() {
      self.terminates++;
      if (self.terminates == self.bars.length) {
        self.terminate();
      }
    }

    return bar;
  },

  terminate: function() {
    this.move(this.bars.length);
    this.stream.clearLine();
    this.stream.cursorTo(0);
  },

  move: function(index) {
    if (!this.stream.isTTY) return;
    this.stream.moveCursor(0, index - this.cursor);
    this.cursor = index;
  },

  tick: function(index, value, options) {
    var bar = this.bars[index];
    if (bar) {
      this.move(index);
      bar.otick(value, options);
    }
  }
};

var mbars = new Multibar();

[quantitative,ambiguous,hexnac_hcd,fragmentation,peptide].forEach(function(module) {
    module.on('task',function(desc) {
        this.addListener('progress',function(percentage) {
            if ( ! progress_bars[module.constructor.name+desc] ) {
                progress_bars[module.constructor.name+desc] = mbars.newBar(module.constructor.name+" "+desc+' [:bar] ', { total: 100} );
            }

            var progress = progress_bars[module.constructor.name+desc];
            progress.tick( (100*percentage).toFixed(0) - progress.curr );
            if (percentage == 1) {
                delete progress_bars[module.constructor.name+desc];
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

var clone = function(objectToBeCloned) {
  // Basis.
  if (!(objectToBeCloned instanceof Object)) {
    return objectToBeCloned;
  }

  var objectClone;

  // Filter out special objects.
  var Constructor = objectToBeCloned.constructor;
  switch (Constructor) {
    // Implement other special objects here.
    case RegExp:
      objectClone = new Constructor(objectToBeCloned);
      break;
    case Date:
      objectClone = new Constructor(objectToBeCloned.getTime());
      break;
    default:
      objectClone = new Constructor();
  }

  // Clone each property.
  for (var prop in objectToBeCloned) {
    objectClone[prop] = clone(objectToBeCloned[prop]);
  }

  return objectClone;
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

var browse_file = function(filename) {
    return open_db(filename).then(function(db) {
        promisify_sqlite(db);
        return db;
    });
};

var process_data = function(filename,sibling_files) {
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
          .then(hexnac_hcd.guess_hexnac.bind(hexnac_hcd,db,sibling_files))
          .then(partial(fragmentation.validate_peptide_coverage,db));

        return processing_promise.then(function(peps) {
            return uniprot_meta.init().then(function() {
                global_datablock.data = peptide.combine(peps);
                quantitative.clear_caches();
                spectra.clear_caches();
                return global_datablock;
            });
        });
    });
};


var files_by_source = {};

var process_files = function(files,sources) {
    var result = Promise.resolve(true);
    var blocks = [];
    files = [].concat(files);
    sources = [].concat(sources);
    files.forEach(function(file,idx) {
        var source = sources[idx];
        files_by_source[source] = (files_by_source[source] || []).concat(file);
    });
    return result.then(function() {
        var self_func = arguments.callee;
        var file = files.shift();
        var source = sources.shift();
        if (! file) {
            return blocks;
        }
        var other_files = files_by_source[source].filter(function(other_file) {
            return other_file !== file;
        });
        return process_data(file,other_files).then(function(block) {
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
                block.data[prot] = block.data[prot].map(function(pep) {
                    return clone(pep);
                });
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

module.exports = exports = { process: process_files, combine: combine, browse_file: browse_file };

