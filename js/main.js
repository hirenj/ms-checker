"use strict";

// Under Node-webkit we need to rebuild the sqlite package

// NODE_WEBKIT_VERSION="0.12.0"
// npm install sqlite3 --build-from-source --runtime=node-webkit --target_arch=x64 --target=$NODE_WEBKIT_VERSION

var gui = null;

try {
    // Load native UI library
    gui = require('nw.gui');
} catch(e) {
    console.log("No GUI available");
}

var nconf = require('nconf');

var fs = require('fs');
var path = require('path');

if (gui) {
    nconf.env({ separator: "__", whitelist : ['MS_DEBUG'] }).add('metadata', { type: 'literal', store: {} }).overrides( require('optimist')(gui.App.argv).argv );

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
} else {
    nconf.env({ separator: "__", whitelist : ['MS_DEBUG'] }).argv();
}

let pretty_print = function(data,depth) {
  let tab_size = 2;
  let indent = Array(tab_size+1).join(' ');
  let output = JSON.stringify(data,null,indent);
  let start_re = new RegExp('\n {'+(depth+1)*tab_size+',}','gm');
  let end_re = new RegExp('([\\"\\}\\]])\\n {'+(depth*tab_size)+'}','gm');
  return output.replace(start_re,'').replace(end_re,'$1');
};



var files_to_open = nconf.get('_') || [];

if (gui) {
    if (files_to_open.length < 1) {
        document.getElementById('fileDialog').click();
    }
}

var sources = nconf.get('source');

if (sources && ! Array.isArray(sources)) {
    sources = [ sources ];
}

var manifest = nconf.get('manifest');

var manifiest_parsing_done = Promise.resolve(true);

if (manifest) {
    files_to_open = [];
    sources = [];

    // Read .tsv or .xls file with
    // data on the manifest contents to build the output files
    manifiest_parsing_done =  Promise.all( [ require('../js/uniprot').cellosaurus.init(), require('../js/entrez').init() ] ).then(function() {
        var conf_data = require('../js/manifests').populateConf(manifest);

        let ki_summary = conf_data['perturbation-ki'].length ? 'KI('+conf_data['perturbation-ki'].map( (ki) => ki.symbol ).join(';')+')' : '';
        let ko_summary = conf_data['perturbation-ko'].length ? 'KO('+conf_data['perturbation-ko'].map( (ko) => ko.symbol ).join(';')+')' : '';
        let other_perturbation_summary = conf_data['perturbation-other'].length > 0 ? 'Other('+conf_data['perturbation-other']+')' : '';
        conf_data.output = [    ( conf_data['source-cell_line'] || (conf_data['source-organism'] + '(' + (conf_data['source-organism-part'] || 'whole') + ')') ),
                                ki_summary, ko_summary, other_perturbation_summary
                            ].filter((val) => val).join('-')+'.msdata';
        nconf.use('metadata', { type: 'literal', store: conf_data });
        files_to_open = nconf.get('input_files').map( (file) =>  path.join(path.dirname(manifest),file) );
        sources = nconf.get('input_sources');
    });
}

var processor;
var excel_writer = require('../js/excel');

var merge_datas = function(current,data) {
    Object.keys(data).forEach(function(up) {
        if ( ! current[up]) {
            current[up] = data[up];
        } else {
            current[up] = current[up].concat(data[up]);
        }
    });
};

// Group files by source
// Check for the ETD parent mass + RT + scan in all the files

manifiest_parsing_done
                     .then(function() { processor = require('../js/processor'); })
                     .then(function() { return processor.process(files_to_open,sources); })
                     .then(function(blocks) { return processor.combine(blocks,sources); })
                     .then(function(combined) {

    // console.log(combined);
    let outpath = nconf.get('output');
    if (nconf.get('outputdir') && outpath) {
        outpath = path.join(nconf.get('outputdir'),nconf.get('output'));
    }
    if (outpath) {
        let existing = '';
        try {
            if ((existing = fs.readFileSync(outpath+'.json')) && nconf.get('manifest')) {
                console.log("Merging existing data for manifest ",nconf.get('manifest'));
                merge_datas(combined.data,JSON.parse(existing).data);
            }
        } catch (err) {
            console.log(err);
        }
        fs.writeFile(outpath+'.json',pretty_print(combined,2),function() {
            console.log("Wrote combined file");
        });
        excel_writer.write(combined,outpath+'.xlsx').catch(console.log.bind(console));
    }
}).catch(function(err) {
    console.error(err);
    console.error(err.stack);
});

