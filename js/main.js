"use strict";

// Under Node-webkit we need to rebuild the sqlite package

// NODE_WEBKIT_VERSION="0.12.0"
// npm install sqlite3 --build-from-source --runtime=node-webkit --target_arch=x64 --target=$NODE_WEBKIT_VERSION

var nconf = require('nconf');


var fs = require('fs');
var path = require('path');

nconf.env({ separator: "__", whitelist : ['MS_DEBUG'] }).argv();

if (nconf.get('help') || ! nconf.get('manifest') || ! nconf.get('_')) {
    console.log("ms-checker --manifest manifestfile.xlsx --outputdir output ..OR.. ms-checker --source somesource path/to/file.msf --output outputfile ");
    process.exit(1);
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
        var path_prefix = '';
        if (nconf.get('prefix-basename')) {
            path_prefix = path.basename(path.dirname(manifest)).replace(/[^A-Za-z0-9]+/,'_').toLowerCase()+'_';
        }

        let wt_summary = conf_data['perturbation-wt'].length ? 'WT('+conf_data['perturbation-wt'].map( (wt) => wt.symbol ).join(';')+')' : '';
        let ki_summary = conf_data['perturbation-ki'].length ? 'KI('+conf_data['perturbation-ki'].map( (ki) => ki.symbol ).join(';')+')' : '';
        let ko_summary = conf_data['perturbation-ko'].length ? 'KO('+conf_data['perturbation-ko'].map( (ko) => ko.symbol ).join(';')+')' : '';
        let other_perturbation_summary = conf_data['perturbation-other'].length > 0 ? 'Other('+conf_data['perturbation-other']+')' : '';
        let source_info = conf_data['source-cell_line'] || conf_data['source-organism'] || '';
        source_info = (''+source_info).replace('-','');
        conf_data.output = path_prefix+[    source_info + '(' + (conf_data['source-organism_part'] || 'bto:0001489') + ')',
                                wt_summary, ki_summary, ko_summary, other_perturbation_summary
                            ].filter((val) => val).join('-')+'.msdata';
        nconf.use('metadata', { type: 'literal', store: conf_data });
        files_to_open = nconf.get('input_files').map( (file) =>  path.join(path.dirname(manifest),file) );
        sources = nconf.get('input_sources');
    });
}

var processor;
var excel_writer = require('../js/excel');
var metadata = require('../js/metadata');

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

var check_output_timestamps = function() {
    var meta = {};
    return metadata.version(meta).then(version => {
        var version_components = meta.software[0].version.split('-');
        var base_version = version_components[0];
        var revision_distance = version_components[1];

        if (nconf.get('use-git-revision')) {
            base_version = base_version + '-' + revision_distance;
        }

        var oldest = new Date(0);

        var files_to_check = [].concat(files_to_open);
        if (nconf.get('manifest')) {
            files_to_check.push(nconf.get('manifest'));
        }

        files_to_check.forEach( file => {
            var file_date = new Date(fs.statSync(file).mtime);
            if (file_date > oldest) {
                oldest = file_date;
            }
        });
        var current_config = nconf.get();
        current_config.output = base_version + '@' + current_config.output;
        nconf.add('metadata',{'type' : 'literal', 'store' : current_config});
        var out_path = path.join(nconf.get('outputdir'),nconf.get('output')+'.json');
        console.log("Deciding if we need to run for ",out_path);
        if (fs.existsSync(out_path)) {
            var current_version_time = (new Date(fs.statSync(out_path).mtime)).getTime();
            if  ( ((new Date()).getTime() - current_version_time) < 5*60*1000 ) {
                console.log("Executing as the existing file is less than 5 minutes old");
                // Run anyway
                return true;
            }
            if ( current_version_time < oldest.getTime() ) {
                console.log("Executing as the target MSF files are newer than the existing file");
                // Run here
                return true;
            }
            console.log("Not running as output files are up-to-date");
            return;
        } else {
            console.log("Running because the target file does not exist");
            // Run here
            return true;
        }
    });
};


// Get the version for the current execution
// Prune version number (conf setting to determine how much pruning)
// Check if file exists in output dir
// If file is < 60 mins old run anyway
// If input MSFS are newer than the output file, run anyway

manifiest_parsing_done
                     .then( check_output_timestamps )
                     .then( (needs_run) => { if ( ! needs_run ) { throw new Error('current') } })
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
            if (err.code !== 'ENOENT') {
                console.log(err);
            }
        }
        fs.writeFile(outpath+'.json',pretty_print(combined,2),function() {
            console.log("Wrote combined file");
        });
        excel_writer.write(combined,outpath+'.xlsx').catch(console.log.bind(console));
    }
}).catch(function(err) {
    if (err.message == 'current') {
        return;
    }
    console.error(nconf.get('input_files'));
    console.error(err);
    console.error(err.stack);
});

