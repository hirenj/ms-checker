// Load native UI library
var gui = require('nw.gui');

var fs      = require('fs'),
    nconf   = require('nconf'),
    sqlite3 = require('sqlite3');

Math.median = require('../js/math-median').median;

var excel_writer = require('../js/excel.js');

var uniprot_meta = require('../js/uniprot.js');

var promisify_sqlite = require('../js/sqlite-promise');

var quantitative = require('../js/quantitative');
var ambiguous = require('../js/ambiguous');
var hexnac_hcd = require('../js/hexnac_hcd');
var fragmentation = require('../js/fragmentions');

var peptide = require('../js/peptide');

ambiguous.conf = nconf;
hexnac_hcd.conf = nconf;

[quantitative,ambiguous,hexnac_hcd,fragmentation].forEach(function(module) {
    module.on('task',function(desc) {
        console.log(module.constructor.name,desc);
        this.addListener('progress',function(percentage) {
            console.log(module.constructor.name,desc, 100*percentage);
            if (percentage == 1) {
                module.removeListener('progress',arguments.callee);
            }
        });
    });
});

var current_files = [];
nconf.env({ separator: "__", whitelist : ['MS_DEBUG'] }).overrides( require('optimist')(gui.App.argv).argv );

if (nconf.get('MS_DEBUG') || nconf.get('debug')) {
    gui.Window.get().showDevTools();
}

(function() {
    var win = gui.Window.get();
    var nativeMenuBar = new gui.Menu({ type: "menubar" });
    try {
        nativeMenuBar.createMacBuiltin("My App");
        win.menu = nativeMenuBar;
    } catch (ex) {
        console.log(ex.message);
    }
})();
var files_to_open = nconf.get('_') || [];
if (files_to_open.length < 1) {
    document.getElementById('fileDialog').click();
}


var partial = function(fn) {
    var args = Array.prototype.slice.call(arguments, 1);
    return function() { return fn.apply(this, args.concat(Array.prototype.slice.call(arguments, 0))); };
};

var onlyUnique = function(value, index, self) {
    return self.indexOf(value) === index;
};

var pbcopy = function(data) { var proc = require('child_process').spawn('pbcopy'); proc.stdin.write(data); proc.stdin.end(); }

var print_pep = function(protein,pep) { return [ protein, pep.sequence, pep.quant ? pep.quant.quant : "", pep.quant ? pep.quant.mad : "" , pep.composition.map(function(comp) { return comp.replace(/x.*/,''); }).join(',') ].join('\t'); };

var filter_global_results = function(datablock,filter) {
    if ( typeof datablock === 'function' ) {
        filter = datablock;
        datablock = global_datablock;
    }
    if ( ! filter ) {
        filter = function() { return true; };
    }
    return "uniprot\tsequence\tquant\tquant_mad\tglyco\n"+
    Object.keys(datablock).map(function(prot) {
        var peps = datablock[prot].filter(filter);
        if (peps.length > 0) {
            return peps.map( partial(print_pep,prot) ).join("\n");
        }
        return "";
    }).filter(function(val) { return val != ''; }).join("\n");
};

var find_protein = function(peps,uniprot) {
    return peps.filter(function(pep) { return pep.uniprot.indexOf(uniprot) >= 0; });
};

var find_sequence = function(peps,sequence) {
    return peps.filter(function(pep) { return pep.Sequence.indexOf(sequence) == 0; });
};

var find_protein_sequence = function(peps,uniprot,sequence) {
    return find_sequence(find_protein(peps,uniprot),sequence);
};


var singlets_filter = function(pep) {
    return pep.quant === 0.00001 || pep.quant === 100000 || pep.quant === 'conflicting_singlets' || pep.quant === 'potential_medium' || pep.quant === 'potential_light';
};

var collect_galnac_ratios = function(peps) {
    if ( ! peps ) {
        peps = global_results;
    }
    return peps.filter(function(pep) { return pep.galnac_intensity > 0; }).map(function(pep) { return pep.galnac_intensity / pep.glcnac_intensity; }).join("\t");
};

var global_results = [];
var global_datablock;

// Load native UI library
var gui = require('nw.gui');

// Print arguments
console.log(gui.App.argv);


var db = new sqlite3.Database(files_to_open[0],sqlite3.OPEN_READONLY,function(err) {

    promisify_sqlite(db);

    var quant_promise = quantitative.init_caches(db).then(function() {
         return quantitative.retrieve_quantified_peptides(db).then(partial(quantitative.check_quantified_peptides,db));
    });
    var ambig_promise = ambiguous.retrieve_peptides(db);

    var processing_promise = Promise.all([quant_promise,ambig_promise]).then(function(all_peps) {
        var merged = Array.prototype.concat.apply([], all_peps);
        return merged;
    }).then(peptide.filter_ambiguous_spectra).then(hexnac_hcd.guess_hexnac.bind(hexnac_hcd,db)).then(partial(fragmentation.validate_peptide_coverage,db));

    processing_promise.then(function(peps) {
        global_results = peps;
        uniprot_meta.init().then(function() {
            console.log(global_datablock);
            global_datablock = peptide.combine(peps);
        });
    }).catch(console.error);

});