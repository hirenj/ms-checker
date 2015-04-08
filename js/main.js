// Load native UI library
var gui = require('nw.gui');

var fs      = require('fs'),
    nconf   = require('nconf'),
    sqlite3 = require('sqlite3'),
    yauzl   = require('yauzl'),
    sax = require('sax'),
    xml2js = require('xml2js'),
    saxpath = require('saxpath');

var heapdump = require('heapdump');


var current_files = [];
nconf.env({ separator: "__", whitelist : ['MS_DEBUG'] }).overrides( require('optimist')(gui.App.argv).argv );

if (nconf.get('MS_DEBUG') || nconf.get('debug')) {
    gui.Window.get().showDevTools();
}

var files_to_open = nconf.get('_') || [];
if (files_to_open.length < 1) {
    document.getElementById('fileDialog').click();
}

var promisify_sqlite = function(db) {
    var old_all = db.all;
    db.all = function(sql,vals) {
        var args = Array.prototype.splice.call(arguments);
        return new Promise(function(resolve,reject) {
            old_all.call(db,sql,vals,function(err,vals) {
                if (err) {
                    reject(err);
                } else {
                    resolve(vals);
                }
            });
        });
    };
    var cached_statements = {};
    db.do_statement = function(sql,vals) {
        if ( ! cached_statements[sql]) {
            cached_statements[sql] = db.prepare(sql);
        }
        return new Promise(function(resolve,reject) {
            cached_statements[sql].all.apply( cached_statements[sql], vals.concat( function(err,vals) {
                if (err) {
                    reject(err);
                } else {
                    resolve(vals);
                }
            }));
        });
    };
    db.end_statement = function(sql) {
        cached_statements[sql].finalize();
        delete cached_statements[sql];
    };
};

const search_peptides_sql = 'SELECT \
    peptides.PeptideID, \
    peptides.UniquePeptideSequenceID, \
    peptides.SpectrumID, \
    peptides.Sequence, \
    PrecursorIonQuanResultsSearchSpectra.QuanResultID, \
    PrecursorIonQuanResults.QuanChannelID, \
    PrecursorIonQuanResults.Area \
FROM peptides \
    LEFT JOIN PrecursorIonQuanResultsSearchSpectra \
        ON peptides.SpectrumID = PrecursorIonQuanResultsSearchSpectra.SearchSpectrumID \
    LEFT JOIN PrecursorIonQuanResults \
        ON PrecursorIonQuanResultsSearchSpectra.QuanResultID = PrecursorIonQuanResults.QuanResultID \
WHERE peptides.ConfidenceLevel = 3 \
AND peptides.PeptideID in (SELECT distinct PeptideID \
    FROM PeptidesAminoAcidModifications \
    WHERE AminoAcidModificationID in (SELECT AminoAcidModificationID \
        FROM AminoAcidModifications \
        WHERE ModificationName like "%Hex%") \
    )';

/*
    We also want to get the subtracted HCD data. We can look at each of the FileNames, and then
    extract out all the MassPeaks (and so Spectra and Peptides) that are associated with each of
    the different filenames, and then assign the ambiguous mass difference to each of the
    peptides.
*/

const related_quants_sql = 'SELECT \
    EventID, \
    FileID, \
    Mass, \
    RT, \
    LeftRT, \
    RightRT, \
    Charge \
FROM EventAnnotations \
    JOIN Events USING(EventID) \
WHERE QuanResultID = ?';

const search_pair_sql = 'SELECT \
    * \
FROM Events \
WHERE   FileID = ? AND \
        RT >= ? AND \
        RT <= ? AND \
        Mass >= ? AND \
        Mass <= ? \
';

const dimethyl_counts_sql = 'SELECT \
    count(distinct Position) as count, PeptideID \
FROM PeptidesAminoAcidModifications \
    JOIN AminoAcidModifications USING(AminoAcidModificationID) \
WHERE (ModificationName = "Dimethyl" OR ModificationName = "Dimethyl:2H(4)") GROUP BY PeptideID';


const peptide_metadata_sql = 'SELECT \
    Description \
FROM Peptides \
    LEFT JOIN PeptidesProteins USING(PeptideID) \
    LEFT JOIN ProteinAnnotations USING (ProteinID) \
WHERE Peptides.PeptideID = ?';

// const peptide_modification_sql = 'SELECT \
//  Position, ModificationName \
// FROM PeptidesAminoAcidModifications \
//  LEFT JOIN AminoAcidModifications USING (AminoAcidModificationID) \
// WHERE PeptideID = ?';

const peptide_modification_sql = 'SELECT \
    Position,AminoAcidModificationID \
FROM PeptidesAminoAcidModifications \
WHERE PeptideID = ?';

const all_peptide_modifications_sql = 'SELECT \
    PeptideID, Position, ModificationName \
FROM PeptidesAminoAcidModifications \
    LEFT JOIN AminoAcidModifications USING (AminoAcidModificationID)';

const retrieve_spectrum_sql = 'SELECT \
    Spectrum \
FROM SpectrumHeaders \
    LEFT JOIN Spectra USING (UniqueSpectrumID) \
WHERE SpectrumID = ? AND SpectrumHeaders.CreatingProcessingNodeNumber = ?';


const retrieve_ambiguous_peptides_sql = 'SELECT \
    Peptides.PeptideID, Peptides.SpectrumID, Peptides.Sequence, FileInfos.Filename \
FROM FileInfos \
    LEFT JOIN MassPeaks USING(FileID) \
    LEFT JOIN SpectrumHeaders USING(MassPeakID) \
    JOIN Peptides using(SpectrumID) \
    WHERE \
        (   FileInfos.FileName LIKE "%xGalNAc%" \
        OR \
            FileInfos.FileName LIKE "%xHexNAc%" \
        OR \
            FileInfos.FileName LIKE "%Hex%" \
        OR \
            FileInfos.FileName LIKE "%Man%" \
        ) \
    AND \
        Peptides.ConfidenceLevel = 3';

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

const MASS_MEDIUM = 32.056407;
const MASS_LIGHT = 28.0313;

var dimethyl_count_cache = null;

var find_dimethyls = function(db,pep) {
    return new Promise(function(resolve,reject) {

        if (! dimethyl_count_cache ) {

            dimethyl_count_cache = {};
            console.log("Populating Dimethyl cache");

            db.all(dimethyl_counts_sql,[]).then(function(data) {
                dimethyl_count_cache = {};
                data.forEach(function(count) {
                    dimethyl_count_cache[count.PeptideID] = count.count+1;
                });
                console.log("Populated Dimethyl cache");

                if ( ! pep.QuanResultID ) {
                    resolve();
                } else {
                    resolve(find_dimethyls(db,pep));
                }
            },reject);

        } else {
            resolve(dimethyl_count_cache[pep.PeptideID]);
            // return check_potential_pair(db,pep,dimethyl_count_cache[pep.PeptideID]);
        }
    });
};


var validated_quans_cache = {};

var check_potential_pair = function(db,pep,num_dimethyl) {
    return new Promise(function(resolve,reject) {
        if (! pep.QuanResultID) {
            resolve();
            return;
        }
        if (pep.QuanResultID in validated_quans_cache) {
            pep.has_pair = validated_quans_cache[pep.QuanResultID];
            resolve();
            return;
        }
        if (pep.Sequence.indexOf('K') < 0 ) {
            num_dimethyl = 1;
        }
        if ( ! num_dimethyl ) {
            find_dimethyls(db,pep).then(function(count) {
                resolve(check_potential_pair(db,pep,count));
            },reject);
            return;
        }

        pep.has_pair = false;

        var mass_change_dir = 1;

        if (pep.QuanChannelID.indexOf(2) >= 0) {
            mass_change_dir = -1;
        }

        db.all(related_quants_sql,[ pep.QuanResultID ]).then(function(events) {
            var events_length = events.length;
            events.forEach(function(ev) {
                events_length -= 1;
                if (pep.has_pair) {
                    if (events_length <= 0) {
                        validated_quans_cache[pep.QuanResultID] = pep.has_pair;
                        resolve();
                    }
                    return;
                }
                var target_mass = ev.Mass + (mass_change_dir * num_dimethyl * (MASS_MEDIUM-MASS_LIGHT)/ev.Charge);
                db.all(search_pair_sql, [ ev.FileID, ev.LeftRT - 0.05, ev.RightRT + 0.05, target_mass*(1 - 15/1000000) , target_mass*(1 + 15/1000000)]).then(function(rows) {
                    if (rows && rows.length > 0) {
                        pep.has_pair = true;
                        if (nconf.get('MS_DEBUG') || nconf.get('debug')) {
                            pep.matching_events = rows.map(function(row) { return row.EventID; });
                            pep.search_event = ev.EventID;
                            pep.target_mass = target_mass;
                            pep.num_dimethyl = num_dimethyl;
                            pep.mass_change_dir = mass_change_dir;
                        }
                    }
                    if (events_length <= 0) {
                        validated_quans_cache[pep.QuanResultID] = pep.has_pair;
                        resolve();
                    }
                }).catch(reject);
            });
        }).catch(reject);
    });
};

var pep_calls = 0;

var produce_peptide_data = function(db,pep) {
    return db.do_statement(peptide_metadata_sql, [ pep.PeptideID ]).then(function(pep_datas) {
        pep_calls++;
        if ((pep_calls % 100) == 0) {
            console.log(pep_calls);
        }
        if ( ! pep.uniprot ) {
            pep.uniprot = [];
        }
        pep_datas.forEach(function(pep_data) {
            var uniprot = pep_data.Description.split('|')[1];
            if (pep.uniprot.indexOf(uniprot) < 0) {
                pep.uniprot.push(uniprot);
            }
        });
    });
};

var peptide_modifications_cache = null;

var produce_peptide_modification_data = function(db,pep) {
    return new Promise(function(resolve,reject) {

        if (! peptide_modifications_cache ) {

            peptide_modifications_cache = {};
            console.log("Populating Modifications cache");

            db.all(all_peptide_modifications_sql,[]).then(function(data) {
                peptide_modifications_cache = {};
                data.forEach(function(mods) {
                    peptide_modifications_cache[mods.PeptideID] =  peptide_modifications_cache[mods.PeptideID] || [];
                    peptide_modifications_cache[mods.PeptideID].push([ mods.Position + 1, mods.ModificationName ]);
                });
                console.log("Populated Modifications cache");

                if ( ! pep.PeptideID ) {
                    resolve();
                } else {
                    resolve(produce_peptide_modification_data(db,pep));
                }
            },reject);

        } else {
            pep.modifications = (peptide_modifications_cache[pep.PeptideID] || []).filter( function(mod) { return (mod[1] || '').indexOf('Hex') >= 0 } ).map( function(mod) { return mod } );
            resolve();
        }
    });
};

var init_caches = function(db) {
    return produce_peptide_modification_data(db,{}).then(function() { return find_dimethyls(db,{'PeptideID': 0}); });
};

var combine_quantified_peptides = function(peps) {
    var all_peps = {};
    peps.forEach(function(pep) {
        var curr_pep = all_peps[pep.PeptideID] ||  pep;
        curr_pep.areas = curr_pep.areas || [];

        if (! Array.isArray(curr_pep.QuanChannelID) ) {
            curr_pep.QuanChannelID =  curr_pep.QuanChannelID ? [ curr_pep.QuanChannelID ] : [];
        }
        if (curr_pep != pep && pep.QuanChannelID) {
            curr_pep.QuanChannelID.push(pep.QuanChannelID);
        }
        if (pep.Area) {
            curr_pep.areas.push(pep.Area);
        }
        all_peps[pep.PeptideID] = curr_pep;
    });
    return Object.keys(all_peps).map(function(key) { return all_peps[key]; });
};

var retrieve_quantified_peptides = function(db) {
    return db.all(search_peptides_sql);
};

var check_quantified_peptides = function(db,peps) {
    console.log("Retrieved peptides to search: "+(peps.length)+" total peptides");

    var combined_peps = combine_quantified_peptides(peps);

    var singlet_peps = combined_peps.filter( function(pep) { return pep.QuanChannelID.filter(onlyUnique).length == 1; } );

    var pair_promises = singlet_peps.map( function(pep) { return check_potential_pair(db,pep,null); } );

    var metadata_promises = combined_peps.map( function(pep) { return produce_peptide_data(db,pep).then( produce_peptide_modification_data(db,pep) );  } );

    return Promise.all(pair_promises.concat(metadata_promises)).then(function() {
        return combined_peps;
    });
};

var validated_hcd_cache = {};

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

var is_galnac_mass = function(mz) {
    var galnac_138 = 138.055;
    var galnac_168 = 168.066;
    return ((mz <= galnac_138*(1 + 15/1000000) && mz >= galnac_138*(1 - 15/1000000)) ||
           (mz <= galnac_168*(1 + 15/1000000) && mz >= galnac_168*(1 - 15/1000000)) );
};

var is_glcnac_mass = function(mz) {
    var glcnac_126 = 126.055;
    var glcnac_144 = 144.065;
    return ((mz <= glcnac_126*(1 + 15/1000000) && mz >= glcnac_126*(1 - 15/1000000)) ||
           (mz <= glcnac_144*(1 + 15/1000000) && mz >= glcnac_144*(1 - 15/1000000)) );
};


var parse_spectrum = function(pep,readStream) {


    if (validated_hcd_cache[pep.SpectrumID]) {
        pep.activation = validated_hcd_cache[pep.SpectrumID].activation;
        if (validated_hcd_cache[pep.SpectrumID].hexnac_type) {
            pep.galnac_intensity = validated_hcd_cache[pep.SpectrumID].galnac_intensity;
            pep.glcnac_intensity = validated_hcd_cache[pep.SpectrumID].glcnac_intensity;
            pep.hexnac_type = validated_hcd_cache[pep.SpectrumID].hexnac_type;
        }
        return true;
    }

    var saxParser = sax.createStream(true);
    var peak_stream = new saxpath.SaXPath(saxParser, '//PeakCentroids/Peak');
    var activation_stream = new saxpath.SaXPath(saxParser, '/MassSpectrum/ScanEvent/ActivationTypes');

    // We should be matching this ActivationType first

    activation_stream.on('match', function(xml) {
        // Check that we're matching a HCD spectrum here
        xml2js.parseString(xml,function(err,result) {
            if (! err ) {
                pep.activation = result.ActivationTypes;
            }
        });
    });

    var galnac_intensity = null;
    var glcnac_intensity = null;

    peak_stream.on('match', function(xml) {
        // Check to see what the ratio (mz-138 + mz-168) / (mz-126 + mz-144) is
        // if it is within 0.4 - 0.6, it is a GalNAc, 2.0 or greater, GlcNAc
        xml2js.parseString(xml,function(err,result) {
            if (! err ) {
                var mass = parseFloat(result.Peak['$'].X);
                var intensity = parseFloat(result.Peak['$'].Y);
                if ( is_glcnac_mass(mass) ) {
                    glcnac_intensity = (glcnac_intensity || 0) + intensity;
                } else if ( is_galnac_mass(mass) ) {
                    galnac_intensity = (galnac_intensity || 0) + intensity;
                }
            }
        });
    });

    var result = new Promise(function(resolve,reject) {
        readStream.on('end', function() {
            validated_hcd_cache[pep.SpectrumID] = { 'activation' : pep.activation };
            if (pep.activation == 'HCD') {
                pep.galnac_intensity = galnac_intensity;
                pep.glcnac_intensity = glcnac_intensity;
                var ratio = galnac_intensity / glcnac_intensity;
                pep.hexnac_type = (ratio <= 0.95) ? 'GalNAc' : ( (ratio >= 1.95) ? 'GlcNAc' : 'Unknown' );
                validated_hcd_cache[pep.SpectrumID].galnac_intensity = galnac_intensity;
                validated_hcd_cache[pep.SpectrumID].glcnac_intensity = glcnac_intensity;
                validated_hcd_cache[pep.SpectrumID].hexnac_type = pep.hexnac_type;
            }
            resolve();
        });
    });

    readStream.pipe(saxParser);

    return result;
};

var init_spectrum_processing_num = function(db) {
    return db.all(hcd_processing_node_number_sql).then(function(nodes) {
        if ( ! nodes ) {
            return;
        }
        var processing = nodes[0].ProcessingNodeNumber;
        return db.all(child_processing_node_numbers_sql,["%,"+processing+",%", "%,"+processing, ""+processing+",%",  processing]).then(function(nodes) {
            if (nodes) {
                return nodes[0].ProcessingNodeNumber;
            }
        });
    });
};

var check_spectrum = function(db,pep,processing_node) {
    return db.all(retrieve_spectrum_sql, [ pep.SpectrumID, processing_node ]).then(function(spectra) {
        if (spectra.length == 0) {
            return;
        }
        var spectrum = spectra[0].Spectrum;
        return unzip_spectrum(spectrum).then(partial(parse_spectrum,pep));
    });
};

var retrieve_ambiguous_peptides = function(db) {
    return db.all(retrieve_ambiguous_peptides_sql).then(function(peps) {
        return Promise.all( peps.map(function(pep) { return produce_peptide_data(db,pep); }) ).then(function() {
            peps.forEach(function(pep) {
                var composition = (pep.FileName || "").match(/\d+\x.+(?=\.mgf)/);
                if ( composition ) {
                    pep.Composition = composition.map(function(comp) { return comp.toString(); });
                }
                delete pep.FileName;
            });
            return peps;
        });
    } );
};


var partial = function(fn) {
    var args = Array.prototype.slice.call(arguments, 1);
    return function() { return fn.apply(this, args.concat(Array.prototype.slice.call(arguments, 0))); };
};

var onlyUnique = function(value, index, self) {
    return self.indexOf(value) === index;
};

var global_results;

var db = new sqlite3.Database(files_to_open[0],sqlite3.OPEN_READONLY,function(err) {

    promisify_sqlite(db);
    init_caches(db).then(function() {
        console.log("Searching Peptides");
        global_db = db;

        // Area is the sum of Areas for the Events asociated with the QuanResultID (from EventAnnotations)
        retrieve_quantified_peptides(db).then(partial(check_quantified_peptides,db)).then(function(peps) {
            db.end_statement(peptide_metadata_sql);
            console.log("Done singlet checking");
            global_results = peps;
        }).then(partial(retrieve_ambiguous_peptides,db)).then(function(peps) {
            console.log("Done ambiguous peptides");
            global_results = global_results.concat(peps);
            return global_results;
        }).then(function(peps) {
            var to_cut = [].concat(peps);
            var processing_node = null;
            var result = init_spectrum_processing_num(db).then(function(node) {
                processing_node = node;
                console.log("We only want spectra from Processing Node ",processing_node);
                return true;
            });
            var split_length = 50;
            result = result.then(function() {
                // heapdump.writeSnapshot();
                console.log("Peps remaining to check : ",to_cut.length);
                if (to_cut.length > 0) {
                    return Promise.all(  to_cut.splice(0,split_length).map(function(pep) { return check_spectrum(db,pep,processing_node); }) ).then(arguments.callee);
                }
                return peps;
            });
            return result;
        }).then(function(peps) {
            console.log("Got all data");
        }).catch(function(err) { console.log(arguments); console.log("Rejected"); });

    });


});