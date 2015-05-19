// Load native UI library
var gui = require('nw.gui');

var nconf = require('nconf');
var processor = require('../js/processor.js');

var excel_writer = require('../js/excel.js');

var fs = require('fs');

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


// Load native UI library
var gui = require('nw.gui');

// Print arguments
console.log(gui.App.argv);

processor.process(files_to_open).then(function(blocks) { return processor.combine(blocks,nconf.get('source')) }).then(function(combined) {
    console.log(combined);
    if (nconf.get('output')) {
        fs.writeFile(nconf.get('output')+'.json',JSON.stringify(combined),function() {
            console.log("Wrote combined file");
        });
        excel_writer.write(combined,nconf.get('output')+'.xls');
    }
}).catch(console.log.bind(console));

