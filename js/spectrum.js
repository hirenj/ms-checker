var yauzl   = require('yauzl'),
    byline  = require('byline'),
    nconf   = require('nconf');

var etd_alias = nconf.get('etd-alias');
var hcd_alias = nconf.get('hcd-alias');

const hcd_processing_node_number_sql = 'SELECT \
    ProcessingNodeNumber \
FROM ProcessingNodeParameters \
WHERE \
    ParameterName = "ActivationTypeFilter" and ParameterValue = "Is#'+(hcd_alias ? hcd_alias : 'HCD')+'"';

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
    Spectrum, Charge, ScanNumbers as scan, RetentionTime as rt, Mass \
FROM SpectrumHeaders \
    LEFT JOIN Spectra USING (UniqueSpectrumID) \
WHERE SpectrumID = ?';


const retrieve_spectrum_from_node_sql = 'SELECT \
    Spectrum, Charge, ScanNumbers as scan, RetentionTime as rt, Mass \
FROM SpectrumHeaders \
    LEFT JOIN Spectra USING (UniqueSpectrumID) \
WHERE SpectrumID = ? AND SpectrumHeaders.CreatingProcessingNodeNumber = ?';


const retrieve_related_spectra_sql = 'SELECT SpectrumID \
FROM SpectrumHeaders \
WHERE MassPeakID IN ( \
    SELECT \
        DISTINCT SpectrumHeaders.MassPeakID \
        FROM SpectrumHeaders \
        WHERE SpectrumID = ? \
)';

const retrieve_related_spectra_by_data_sql = 'SELECT SpectrumID,ScanNumbers \
FROM SpectrumHeaders \
WHERE \
    RetentionTime BETWEEN ? AND ? \
    AND Charge = ? \
    AND Mass BETWEEN ? AND ?';

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

var fix_activation = function(spectrum) {
    if (etd_alias || hcd_alias) {
        if (etd_alias && spectrum.activation == etd_alias) {
            spectrum.activation = 'ETD';
            return spectrum;
        }
        if (hcd_alias && spectrum.activation == hcd_alias) {
            spectrum.activation = 'HCD';
            return spectrum;
        }
        return spectrum;
    }

    // Retrieve the mapping of ETD and HCD identifier
    // from the params and make this function a no-op if
    // we don't have any mapping we need to do

    if (nconf.get('hcd-alias')) {
        hcd_alias = nconf.get('hcd-alias');
    }
    if (nconf.get('etd-alias')) {
        etd_alias = nconf.get('etd-alias');
    }

    if ( ! hcd_alias && ! etd_alias ) {
        fix_activation = function(spec) { return spec };
    }

    return fix_activation(spectrum);
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
            if (line.match('<InstrumentName>')) {
                spectrum_object.instrument = line.split(/[<>]/)[2];
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
    var instrument_stream = new saxpath.SaXPath(saxParser, '/MassSpectrum/Header/InstrumentName');

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

    instrument_stream.on('match',function(xml) {
        xml2js.parseString(xml,function(err,result) {
            if (! err ) {
                spectrum_object.instrument = result.InstrumentName;
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

var get_spectrum = function(db,pep,processing_node,no_cache) {
    var spectrum_promise;
    if (processing_node === 0) {
        processing_node = null;
    }
    var cache_id = pep.SpectrumID;
    if (no_cache) {
        cache_id = no_cache+pep.SpectrumID;//Math.random().toString(36).substring(10);
    }

    if ( ! spectrum_caches[cache_id] ) {
        if ( typeof processing_node !== 'undefined' && processing_node !== null ) {
            spectrum_promise = db.all(retrieve_spectrum_from_node_sql, [ pep.SpectrumID, processing_node ]);
        } else {
            spectrum_promise = db.all(retrieve_spectrum_sql, [ pep.SpectrumID ]);
        }
    }
    spectrum_caches[cache_id] = spectrum_caches[cache_id] || spectrum_promise.then(function(spectra) {
        if (spectra.length == 0) {
            delete spectrum_caches[cache_id];
            return;
        }
        var spectrum = spectra[0].Spectrum;
        var charge = spectra[0].Charge;
        var mass = spectra[0].Mass;
        var scan = spectra[0].scan;
        var retentionTime = spectra[0].rt;
        return unzip_spectrum(spectrum).then(parse_spectrum).then(fix_activation).then(function(spectrum){
            spectrum.spectrumID = pep.SpectrumID;
            spectrum.charge = charge;
            spectrum.mass = mass;
            spectrum.scan = scan;
            spectrum.rt = retentionTime;
            if ( ! db.instrument && spectrum.instrument ) {
                db.instrument = spectrum.instrument;
            }
            if ( ! spectrum.instrument ) {
                spectrum.instrument = db.instrument;
            }
            return spectrum;
        });
    });

    spectrum_caches[cache_id].then(function(spectrum) {
        pep.retentionTime = spectrum.rt;
        pep.scan = spectrum.scan;
        pep.mass = spectrum.mass;
        pep.charge = spectrum.charge;
        pep.activation = spectrum.activation;
    });

    return spectrum_caches[cache_id];
};

var get_related_spectra = function(db,spectrum) {
    if ( ! spectrum.spectrumID ) {
        return [];
    }
    return db.all(retrieve_related_spectra_sql,[spectrum.spectrumID]).then(function(spec_ids) {
        return Promise.all((spec_ids || []).map(function(spec_id) {
            return get_spectrum(db,spec_id);
        }));
    });
};

var filter_hcd_with_mods = function(peps) {
    peps = peps.filter(function(pep) {
        return (pep.activation !== 'HCD' && pep.activation !== 'CID') || ! pep.modifications;
    });
    return peps;
};

var match_spectrum_data = function(db, scan, rt, charge, mass) {
    var wanted_scan = parseInt(scan);
    return db.do_statement(retrieve_related_spectra_by_data_sql,[rt - 0.03,rt + 0.03,charge,mass-0.01,mass+0.01]).then(function(spec_ids) {
        var wanted_spec_ids = (spec_ids || []).filter(function(spec) {  return Math.abs(wanted_scan - parseInt(spec.ScanNumbers)) <= 3;  });
        return Promise.all(wanted_spec_ids.map(function(spec_id) {
            return get_spectrum(db,spec_id,null,db.partner);
        }));
    });
};

var clear_caches = function() {
    spectrum_caches = {};
};



exports.init_spectrum_processing_num = init_spectrum_processing_num;
exports.get_spectrum = get_spectrum;
exports.get_related_spectra = get_related_spectra;
exports.match_spectrum_data = match_spectrum_data;
exports.clear_caches = clear_caches;
exports.filter_hcd_with_mods = filter_hcd_with_mods;

