
var nconf = require('nconf');
var sqlite3 = require('sqlite3');
var yauzl   = require('yauzl'),
    byline  = require('byline');

var fs = require('fs');
nconf.env({ separator: "__", whitelist : ['MS_DEBUG'] }).argv();

var promisify_sqlite = require('../js/sqlite-promise');
var spectra = require('../js/spectrum');

var open_db = function(filename) {
    return new Promise(function(resolve,reject) {
        var db = new sqlite3.Database(filename,sqlite3.OPEN_READWRITE,function(err) {
            if (err) {
                reject(err);
                return;
            }
            promisify_sqlite(db);
            resolve(db);
        });
    });
};

var unzip_spectrum = function(spectrum_text) {
    return new Promise(function(resolve,reject) {
        yauzl.fromBuffer(spectrum_text, function(err, zipfile) {

            if (err) {
                reject(err);
                return;
            }
            zipfile.addListener('end',resolve);
            zipfile.addListener("entry", function(entry) {
                if (/\/$/.test(entry.fileName)) {
                    // directory file names end with '/'
                    return;
                }
                zipfile.removeListener('entry',arguments.callee);
                zipfile.removeListener('end',resolve);

                zipfile.openReadStream(entry, function(err, readStream) {
                    if (err)  {
                        reject(err);
                        return;
                    }
                    resolve(readStream);
                });
            });
        });
    });

};


const retrieve_spectrum_sql = 'SELECT \
    Spectrum, Charge, ScanNumbers as scan, RetentionTime as rt, Mass \
FROM SpectrumHeaders \
    LEFT JOIN Spectra USING (UniqueSpectrumID) \
ORDER BY RANDOM() LIMIT 1';

const retrieve_spectrum_sql_fixed = 'SELECT \
    Spectrum, Charge, ScanNumbers as scan, RetentionTime as rt, Mass \
FROM SpectrumHeaders \
    LEFT JOIN Spectra USING (UniqueSpectrumID) \
WHERE SpectrumHeaders.SpectrumID = ?';

const spec_ids = nconf.get('_').length > 1 ? [ nconf.get('_')[1] ] : [];

const spectrum_sql = spec_ids.length > 0 ? retrieve_spectrum_sql_fixed : retrieve_spectrum_sql;

open_db(nconf.get('_')[0]).then(function(db) {
  db.all(spectrum_sql, spec_ids).then(function(spectra) {
    console.log(spectra[0].Spectrum);
    unzip_spectrum(spectra[0].Spectrum).then(function(vals) {
      vals.on('data',function(data) {
        console.log(data.toString());
      });
    }).catch(function(err) {
      console.error(err);
    });
  });
});