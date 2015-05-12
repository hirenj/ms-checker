var yauzl   = require('yauzl'),
    byline  = require('byline');

const hcd_processing_node_number_sql = 'SELECT \
    ProcessingNodeNumber \
FROM ProcessingNodeParameters \
WHERE \
    ParameterName = "ActivationTypeFilter" and ParameterValue = "Is#HCD"';

const child_processing_node_numbers_sql = 'SELECT ProcessingNodeNumber \
FROM ProcessingNodes \
WHERE (ProcessingNodeParentNumber LIKE ? \
    OR ProcessingNodeParentNumber LIKE ? \
    OR ProcessingNodeParentNumber LIKE ? \
    OR ProcessingNodeParentNumber = ?) \
    AND ProcessingNodeNumber IN ( \
        SELECT \
            distinct CreatingProcessingNodeNumber \
        FROM SpectrumHeaders \
        )';

const retrieve_spectrum_sql = 'SELECT \
    Spectrum, Charge \
FROM SpectrumHeaders \
    LEFT JOIN Spectra USING (UniqueSpectrumID) \
WHERE SpectrumID = ?';


const retrieve_spectrum_from_node_sql = 'SELECT \
    Spectrum, Charge \
FROM SpectrumHeaders \
    LEFT JOIN Spectra USING (UniqueSpectrumID) \
WHERE SpectrumID = ? AND SpectrumHeaders.CreatingProcessingNodeNumber = ?';


var spectrum_caches = {};

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

var parse_spectrum = function(readStream) {
    var stream = byline.createStream(readStream);

    var in_scan_event = false;
    var in_peak_centroids = false;
    var spectrum_object = {};
    spectrum_object.peaks = [];

    var result = new Promise(function(resolve,reject) {
        readStream.on('end', function() {
            resolve( spectrum_object );
        });
    });

    stream.on('readable', function() {
        var line;
        while (null !== (line = stream.read())) {
            line = line.toString().trim();
            if (line == '<ScanEvent>') {
                in_scan_event = true;
            }
            if (line == '</ScanEvent>') {
                in_scan_event = false;
            }
            if (! spectrum_object.activation && in_scan_event && line.indexOf('ActivationTypes') >= 0 ) {
                spectrum_object.activation = line.match(/>(.+)</)[1];
            }
            if (line == '<PeakCentroids>') {
                in_peak_centroids = true;
            }
            if (line == '</PeakCentroids>') {
                in_peak_centroids = false;
            }
            if (in_peak_centroids && line.indexOf('Peak ') >= 0) {
                var peak = {};
                line.split(/\s+/).filter(function(attr) { return attr.indexOf('=') >= 0; }).forEach(function(attr) {
                    var vals = attr.split('=');
                    var attname = vals[0];
                    var value = vals[1];
                    value = value.replace(/"/g,'');
                    if (attname == 'X') {
                        peak.mass = parseFloat(value);
                    }
                    if (attname == 'Y') {
                        peak.intensity = parseFloat(value);
                    }
                    if (attname == 'Z') {
                        peak.charge = parseFloat(value);
                    }
                    if (attname == 'SN') {
                        peak.sn = parseFloat(value);
                    }
                });
                spectrum_object.peaks.push(peak);
            }
        }
    });

    return result;
};


var parse_spectrum_xml = function(readStream) {

    var sax     = require('sax'),
    xml2js  = require('xml2js'),
    saxpath = require('saxpath');

    var saxParser = sax.createStream(true);
    var peak_stream = new saxpath.SaXPath(saxParser, '//PeakCentroids/Peak');
    var activation_stream = new saxpath.SaXPath(saxParser, '/MassSpectrum/ScanEvent/ActivationTypes');

    // ActivationType normally is the first thing matched

    var spectrum_object = {};

    activation_stream.on('match', function(xml) {
        // Check that we're matching a HCD spectrum here
        xml2js.parseString(xml,function(err,result) {
            if (! err ) {
                spectrum_object.activation = result.ActivationTypes;
            }
        });
    });

    spectrum_object.peaks = [];

    peak_stream.on('match', function(xml) {
        xml2js.parseString(xml,function(err,result) {
            if (! err ) {
                spectrum_object.peaks.push({
                    mass: parseFloat(result.Peak['$'].X),
                    intensity: parseFloat(result.Peak['$'].Y),
                    charge: parseFloat(result.Peak['$'].Z),
                    sn: parseFloat(result.Peak['$'].SN)
                });
            }
        });
    });

    var result = new Promise(function(resolve,reject) {
        readStream.on('end', function() {
            resolve( spectrum_object );
        });
    });

    readStream.pipe(saxParser);

    return result;
};

var init_spectrum_processing_num = function(db) {
    return db.all(hcd_processing_node_number_sql).then(function(nodes) {
        if ( ! nodes || nodes.length < 1 ) {
            return;
        }
        var processing = nodes[0].ProcessingNodeNumber;
        return db.all(child_processing_node_numbers_sql,["%,"+processing+",%", "%,"+processing, ""+processing+",%",  processing]).then(function(nodes) {
            if (nodes && nodes.length > 0) {
                return nodes[0].ProcessingNodeNumber;
            }
        });
    });
};

var get_spectrum = function(db,pep,processing_node) {
    var spectrum_promise;
    if ( ! spectrum_caches[pep.SpectrumID] ) {
        if ( typeof processing_node !== 'undefined' && processing_node !== null ) {
            spectrum_promise = db.all(retrieve_spectrum_from_node_sql, [ pep.SpectrumID, processing_node ]);
        } else {
            spectrum_promise = db.all(retrieve_spectrum_sql, [ pep.SpectrumID ]);
        }
    }
    spectrum_caches[pep.SpectrumID] = spectrum_caches[pep.SpectrumID] || spectrum_promise.then(function(spectra) {
        if (spectra.length == 0) {
            delete spectrum_caches[pep.SpectrumID];
            return;
        }
        var spectrum = spectra[0].Spectrum;
        var charge = spectra[0].Charge;
        return unzip_spectrum(spectrum).then(parse_spectrum).then(function(spectrum){
            spectrum.spectrumID = pep.SpectrumID;
            spectrum.charge = charge;
            return spectrum;
        });
    });
    return spectrum_caches[pep.SpectrumID];
};

exports.init_spectrum_processing_num = init_spectrum_processing_num;
exports.get_spectrum = get_spectrum;