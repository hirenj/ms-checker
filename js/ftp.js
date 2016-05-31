var ftp = require('ftp');
var parse_url = require('url').parse;
var fs = require('fs');

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
                resolve([ filename, true ]);
            });
            stream.pipe(fs.createWriteStream(filename));
        });
    });
};

var check_modified = function(timedata,url,filename,grace_period) {
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
            if ( ((new Date(fsstat.mtime)).getTime()+grace_period) < date.getTime() ) {
                resolve( download_file( ftp_site, parse_url(url).path, filename ).then(function(filename) {
                    fs.utimesSync(filename, date,date);
                }) );
            } else {
                ftp_site.end();
                resolve( [filename, false] );
            }
        });
    });
    return result;
};

var get_cached_file = function(url,filename,grace_period) {
    return Promise.all( [ get_modtime(filename), connect_ftp(url) ] ).then(function(promise_results) {
        return check_modified(promise_results,url,filename,grace_period || 0);
    });
};

exports.get_cached_file = get_cached_file;
