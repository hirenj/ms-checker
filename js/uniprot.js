var util = require('util');
var ftp = require('ftp');
var parse_url = require('url').parse;
var fs = require('fs');

// Uniprot deleted accessions (find deleted entries):
const DELETED_URL = 'ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/docs/delac_sp.txt';

// Uniprot secondary accessions (find merged entries):
const MERGED_URL = 'ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/docs/sec_ac.txt';

// Cellosaurus respository
const CELLOSAURUS_URL = 'ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt';

var UniprotMeta = function UniprotMeta() {

};

util.inherits(UniprotMeta,require('events').EventEmitter);

var CellosaurusMeta = function CellosaurusMeta() {

};

util.inherits(CellosaurusMeta,require('events').EventEmitter);

exports = module.exports = {'uniprot': new UniprotMeta(), 'cellosaurus' : new CellosaurusMeta() };

var get_modtime = function(filename) {
    return new Promise(function(resolve,reject) {
        try {
            fs.stat(filename,function(err,stats) {
                if (err) {
                    reject(err);
                    return;
                }
                resolve(stats);
            });
        } catch(e) {
            reject(e);
        }
    }).catch(function(err) { return Promise.resolve({'mtime' : '1-1-1'} ); });
};

var connect_ftp = function(url) {
    var req = parse_url(url);
    var ftp_site  = new ftp();
    var result = new Promise(function(resolve,reject) {
        ftp_site.on('ready',function() {
            resolve(ftp_site);
        });
        ftp_site.on('error',reject);
    });
    ftp_site.connect(req);
    return result.catch(function(err) {
        console.log(err);
        return null;
    });
};

var download_file = function(ftp,path,filename) {
    return new Promise(function(resolve,reject) {
        ftp.get(path, function(err, stream) {
            if (err) throw err;
            stream.once('close', function() {
                ftp.end();
                resolve(filename);
            });
            stream.pipe(fs.createWriteStream(filename));
        });
    });
};

var check_modified = function(timedata,url,filename) {
    var fsstat = timedata[0];
    var ftp_site = timedata[1];

    if ( ! ftp_site ) {
        return Promise.resolve(filename);
    }

    var result = new Promise(function(resolve,reject) {
        ftp_site.lastMod(parse_url(url).path,function(err,date) {
            if (err) {
                ftp_site.end();
                reject(err);
            }
            if ((new Date(fsstat.mtime)).getTime() < date.getTime() ) {
                resolve( download_file( ftp_site, parse_url(url).path, filename ).then(function(filename) {
                    fs.utimesSync(filename, date,date);
                }) );
            } else {
                ftp_site.end();
                resolve( filename );
            }
        });
    });
    return result;
};

var get_cached_file = function(url,filename) {
    return Promise.all( [ get_modtime(filename), connect_ftp(url) ] ).then(function(promise_results) {
        return check_modified(promise_results,url,filename);
    });
};

var metadata = null;

var parse_metadata = function(filename) {
    var uniprot_re = /[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}/;
    return new Promise(function(resolve,reject) {
        fs.readFile(filename, function(err,data){
            console.log(filename);
            if (err){
                reject(err);
            }
            data.toString().split('\n').filter(function(line) { return line.match(uniprot_re); }).forEach(function(line) {
                var bits = line.split(/\s+/);
                metadata[bits[0]] = bits[1] ? bits[1] : -1;
            });
            resolve();
        });
    });
};

var parse_cellosaurus = function(filename) {
    return new Promise(function(resolve,reject) {
        var current_entry = { names : [] };
        cellosaurus = {};
        fs.readFile(filename, function(err,data){
            if (err){
                reject(err);
            }
            data.toString().split('\n').filter(function(line) { return line.match(/(^ID|^AC|^SY|^OX|^DR|^HI|^\/\/)/); }).forEach(function(line) {
                if (line == '//') {
                    if (current_entry) {
                        current_entry.names.forEach(function(name) {
                            cellosaurus[name] = current_entry;
                        });
                    }
                    current_entry = { names: [] };
                    return;
                }
                var line_id = line.substr(0,5).trim();
                var line_data = line.substr(5).trim();
                if (line_id == 'ID' || line_id == 'SY') {
                    current_entry.names = current_entry.names.concat(line_data.split(/;\s+/));
                }
                if (line_id == 'AC') {
                    current_entry.acc = line_data;
                }
                if (line_id == 'OX') {
                    current_entry.taxid = parseInt((/NCBI_TaxID=(\d+)/.exec(line_data) || [undefined,undefined])[1]);
                }
                if (line_id == 'DR' && line_data.match(/BTO\;/)) {
                    current_entry.bto = (/BTO:(\d+)/.exec(line_data) || [undefined,undefined])[1];
                }
                if (line_id == 'HI') {
                    current_entry.parent = (/!(.*)/.exec(line_data) || [undefined,undefined])[1].trim();
                }
            });
            Object.keys(cellosaurus).forEach(function(name) {
                var entry = cellosaurus[name];
                if (entry.parent && ! entry.bto) {
                    entry.bto = cellosaurus[entry.parent].bto;
                }
            });
            console.log(cellosaurus['CHO-GS']);
            resolve(true);
        });
    });
};

UniprotMeta.prototype.init = function() {
    if (metadata) {
        return Promise.resolve(true);
    }
    return Promise.all( [   get_cached_file( DELETED_URL, 'uniprot.deleted.txt' ),
                            get_cached_file( MERGED_URL , 'uniprot.merged.txt' ) ]).then(function(filenames) {
        metadata = {};
        return Promise.all( filenames.map(parse_metadata) ).catch(function(err) { metadata = null; throw err; });
    });
};

UniprotMeta.prototype.isReplaced = function(uniprot) {
    if (! metadata ) {
        return;
    }
    return metadata[uniprot];
};

UniprotMeta.prototype.updateIds = function(ids) {
    var self = this;
    return (ids || []).map(function(id) {
        return self.isReplaced(id) || id;
    }).filter(function(id) { return id !== -1; });
};

var cellosaurus = null;

CellosaurusMeta.prototype.init = function() {
    console.log("Initing cellosaurus");
    if (cellosaurus) {
        return Promise.resolve(true);
    }
    return get_cached_file( CELLOSAURUS_URL, 'cellosaurus.txt').then(function(filename) {
        cellosaurus = {};
        return parse_cellosaurus(filename).catch(function(err) { cellosaurus = null; throw err; });
    });
};

CellosaurusMeta.prototype.lookup = function(name) {
    return Object.freeze(cellosaurus[name]);
};
