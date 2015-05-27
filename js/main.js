var gui = null;

try {
    // Load native UI library
    gui = require('nw.gui');
} catch(e) {
    console.log("No GUI available");
}

var nconf = require('nconf');
var processor = require('../js/processor');

var excel_writer = require('../js/excel');

var fs = require('fs');

if (gui) {
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
}
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

// Group files by source
// Check for the ETD parent mass + RT + scan in all the files

processor.process(files_to_open,sources).then(function(blocks) { return processor.combine(blocks,sources) }).then(function(combined) {
    console.log(combined);
    if (nconf.get('output')) {
        fs.writeFile(nconf.get('output')+'.json',JSON.stringify(combined),function() {
            console.log("Wrote combined file");
        });
        excel_writer.write(combined,nconf.get('output')+'.xls').catch(console.log.bind(console));
    }
}).catch(console.log.bind(console));

