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
        if (cached_statements[sql]) {
            cached_statements[sql].finalize();
        }
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
AND SearchEngineRank = 1 \
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
    SN, \
    QuanChannelID, \
    Charge \
FROM EventAnnotations \
    JOIN Events USING(EventID) \
WHERE QuanResultID = ? AND QuanChannelID = ?';

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
    PeptideID, Position, ModificationName, DeltaMass \
FROM PeptidesAminoAcidModifications \
    LEFT JOIN AminoAcidModifications USING (AminoAcidModificationID) \
UNION \
\
SELECT \
    PeptideID, 0, ModificationName, DeltaMass \
FROM PeptidesTerminalModifications \
    LEFT JOIN AminoAcidModifications ON PeptidesTerminalModifications.TerminalModificationID = AminoAcidModifications.AminoAcidModificationID \
';



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

const MASS_H = 1.007825;
const MASS_C = 12.000000;
const MASS_N = 14.003070;
const MASS_O = 15.994910;

const MASS_AMINO_ACIDS = {
    A:71.03712,
    C:103.00919,
    D:115.02695,
    E:129.0426,
    F:147.06842,
    G:57.02147,
    H:137.05891,
    I:113.08407,
    K:128.09497,
    L:113.08407,
    M:131.0405,
    N:114.04293,
    P:97.05277,
    Q:128.05858,
    R:156.10112,
    S:87.03203,
    T:101.04768,
    V:99.06842,
    W:186.07932,
    Y:163.06332
};


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
        var sns = [];
        pep.has_low_sn = false;

        db.all(related_quants_sql,[ pep.QuanResultID, pep.QuanChannelID.filter(onlyUnique)[0] ]).then(function(events) {
            var events_length = events.length;
            events.forEach(function(ev) {

                sns.push(ev.SN);

                events_length -= 1;
                if (events_length <= 0) {
                    if (median(sns) < 10) {
                        pep.has_low_sn = true;
                    }
                }
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
            if (uniprot && pep.uniprot.indexOf(uniprot) < 0) {
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
                    peptide_modifications_cache[mods.PeptideID].push([ mods.Position + 1, mods.ModificationName, mods.DeltaMass ]);
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
    spectrum_caches[pep.SpectrumID] = spectrum_caches[pep.SpectrumID] || db.all(retrieve_spectrum_sql, [ pep.SpectrumID, processing_node ]).then(function(spectra) {
        if (spectra.length == 0) {
            return;
        }
        var spectrum = spectra[0].Spectrum;
        return unzip_spectrum(spectrum).then(parse_spectrum).then(function(spectrum){
            spectrum.spectrumID = pep.spectrumID;
            return spectrum;
        });
    });
    return spectrum_caches[pep.SpectrumID];
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

var check_galnac_glcnac_ratio = function(pep,spectrum) {
    if (! spectrum ) {
        return;
    }
    var galnac_intensity = null;
    var glcnac_intensity = null;
    var galnac_count = 0;
    var glcnac_count = 0;
    if (spectrum.activation !== 'HCD') {
        return;
    }
    pep.activation = 'HCD';

    // Check to see what the ratio (mz-138 + mz-168) / (mz-126 + mz-144) is
    // if it is within 0.4 - 0.6, it is a GalNAc, 2.0 or greater, GlcNAc

    spectrum.peaks.forEach(function(peak) {
        var mass = peak.mass;
        var intensity = peak.intensity;
        if ( is_glcnac_mass(mass) ) {
            glcnac_intensity = (glcnac_intensity || 0) + intensity;
            glcnac_count += 1;
        } else if ( is_galnac_mass(mass) ) {
            galnac_intensity = (galnac_intensity || 0) + intensity;
            galnac_count += 1;
        }
    });

    if (galnac_count > 2 || glcnac_count > 2) {
        console.log("More peaks counted than wanted for ",pep.SpectrumID);
    }

    if (galnac_count == 2 && glcnac_count == 2) {
        pep.galnac_intensity = galnac_intensity;
        pep.glcnac_intensity = glcnac_intensity;
        var ratio = galnac_intensity / glcnac_intensity;
        pep.hexnac_type = (ratio <= 0.95) ? 'GalNAc' : ( (ratio >= 1.95) ? 'GlcNAc' : 'Unknown' );
    }

    return;
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
//https://github.com/compomics/thermo-msf-parser/blob/master/thermo_msf_parser_API/src/main/java/com/compomics/thermo_msf_parser_API/highmeminstance/Peptide.java
var assign_peptide_ions = function(pep,spectrum) {

};

// calculate_fragment_ions(global_results.filter(function(pep) { return pep.Sequence == 'YSEFFTGSK'; })[0],1);

var calculate_fragment_ions = function(pep,spectrum_charge) {
    var mods = peptide_modifications_cache[pep.PeptideID];
    var null_mass = [0,0,0];
    var masses = pep.Sequence.split('').map(function(aa,idx) { return MASS_AMINO_ACIDS[aa] + (mods.filter(function(mod) { return mod[0] === (idx + 1); })[0] || null_mass)[2];  });
    var ions = [];
    var charge = spectrum_charge;
    var b_ion_base, y_ion_base;
    while(charge > 0) {
        b_ion_base = 0;
        y_ion_base = 2*MASS_H + MASS_O;

        for (var i = 0 ; i < pep.Sequence.length; i++) {
            b_ion_base += masses[i];
            y_ion_base += masses[(masses.length - 1) - i];

            // FIXME: Check calculation for b_nh3 ion
            // { 'type' : 'b_nh3_'+(i+1), 'mz' : (b_ion_base - MASS_O + 2 * MASS_N - 3* MASS_H + charge * MASS_H) / charge, 'z' : charge },

            ions = ions.concat( [
            { 'type' : 'b'+(i+1), 'mz' : (b_ion_base + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'b_h2o_'+(i+1), 'mz' : (b_ion_base - MASS_O - 2 * MASS_H + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'a'+(i+1), 'mz' : (b_ion_base - MASS_O - MASS_C + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'a_nh3_'+(i+1), 'mz' : (b_ion_base - MASS_O - MASS_C - MASS_N - 3 * MASS_H + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'a_h20_'+(i+1), 'mz' : (b_ion_base - 2*MASS_O - MASS_C - 2 * MASS_H + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'c'+(i+1), 'mz' : (b_ion_base + MASS_N + 3 * MASS_H + charge * MASS_H) / charge, 'z' : charge },

            { 'type' : 'y'+(i+1), 'mz' : (y_ion_base + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'y_nh3_'+(i+1), 'mz' : (y_ion_base - MASS_N - 3 * MASS_H + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'y_h2o_'+(i+1), 'mz' : (y_ion_base - 2* MASS_H - MASS_O + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'x'+(i+1), 'mz' : (y_ion_base + MASS_C + MASS_O - 2 * MASS_H + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'z'+(i+1), 'mz' : (y_ion_base - MASS_N - 2 * MASS_H + charge * MASS_H) / charge, 'z' : charge }
            ]);
        }
        charge -= 1;
    }
    return ions;
};

var modification_key = function(pep) {
    var mod_key = null;
    if (pep.modifications) {
        mod_key = pep.modifications.sort(function(a,b) { return a[0] - b[0]; }).map(function(mod) {  return mod[0]+"-"+mod[1]; }).join(',');
    }
    if (pep.Composition) {
        mod_key = pep.Composition.join(',');
    }
    return mod_key;
};

var median = function(values) {

    values.sort( function(a,b) {return a - b;} );

    var half = Math.floor(values.length/2);

    if(values.length % 2)
        return values[half];
    else
        return (values[half-1] + values[half]) / 2.0;
}

var not_hcd_filter = function(pep) {
    return pep.activation !== 'HCD';
};

var combine_all_peptides = function(peps) {
    var data = {};
    var peptides_by_spectrum = {};
    peps.forEach(function(pep) {
        if ( ! peptides_by_spectrum[pep.SpectrumID]) {
            peptides_by_spectrum[pep.SpectrumID] = [];
        }
        peptides_by_spectrum[pep.SpectrumID].push(pep);
    });

    // Search by Spectrum ID:
    // Mark ambiguous peptides - if we have different sets of modifications for each spectrum, we can count that as ambiguous. Maybe mark out the possible amino acids? + total number of sites
    Object.keys(peptides_by_spectrum).forEach(function(spectrumID) {

        // Don't accept any peptides that has been identified using HCD
        // we don't trust the assignment, so we can skip them

        var peps = peptides_by_spectrum[spectrumID];

        // We want to count the number of modification keys for Non-HCD spectra
        // that identify sites. If there are multiple modification keys
        // it's a non-unique identification, so we should just use the
        // composition

        var mods = peps.filter(not_hcd_filter).map(function(pep) { return modification_key(pep); }).filter(onlyUnique).filter(function(key) { return key.indexOf('-') >= 0; });

        peps.forEach(function(pep) {
            if (mods.length > 1 || pep.activation == 'HCD') {
                if ( ! pep.Composition && pep.modifications.length > 0 ) {
                    // We should be a bit smarter about this, grouping compositions appropriately
                    pep.Composition =  [ ""+pep.modifications.length + "x" + pep.modifications[0][1] ];
                }
                pep.possible_mods = pep.modifications;
                delete pep.modifications;
            }
        });

    });
    peptides_by_spectrum = null;

    var grouped_peptides = {};
    peps.forEach(function(pep) {
        var pep_key = pep.Sequence+"-"+modification_key(pep);
        if ( ! grouped_peptides[pep_key] ) {
            grouped_peptides[pep_key] = [];
        }
        grouped_peptides[pep_key].push(pep);
    });

    // Search by Peptide / (modifications + ambig key)
    // Combine quantitated peptides - seperating out into ambiguous / unambiguous.
    Object.keys(grouped_peptides).forEach(function(pep_key) {

        // HCD identification is not used for quantification

        var peps = grouped_peptides[pep_key].filter(not_hcd_filter);
        var quan_peps = peps.filter(function(pep) { return "QuanResultID" in pep;  });
        var singlet_peps = quan_peps.filter(function(pep) { return pep.QuanChannelID.length < 2 && pep.has_pair == false; });
        var ratioed_peps = quan_peps.filter(function(pep) { return pep.QuanChannelID.length == 2; });
        var target_ratio = null;
        if (ratioed_peps.length > 0) {
            var seen_quan_result_ids = {};

            // We wish to consider each ratio once for each
            // quantitation result. If there are multiple modifications
            // or some other thing that boosts the number of peptides
            // this will ignore them.

            target_ratio = median( ratioed_peps.filter(function(pep) {
                var seen = pep.QuanResultID in seen_quan_result_ids;
                seen_quan_result_ids[pep.QuanResultID] = 1;
                return ! seen;
            }).map(function(pep) { return pep.areas[1] / pep.areas[0]; }) );

            singlet_peps.forEach(function(pep) { pep.used = false; });
            ratioed_peps.forEach(function(pep) { pep.used = true; });

        } else if (singlet_peps.length > 0)  {

            target_ratio = singlet_peps[0].QuanChannelID[0] == 1 ? 1/100000 : 100000;

            var channel_ids = singlet_peps.map(function(pep) { return pep.QuanChannelID[0]; }).filter(onlyUnique);

            if (channel_ids.length > 1) {
                console.log(singlet_peps);
                target_ratio = 'conflicting_singlets';
            }
        }
        if (target_ratio) {
            quan_peps.forEach(function(pep) { pep.CalculatedRatio = target_ratio; });
        }
    });

    Object.keys(grouped_peptides).forEach(function(pep_key) {
        var peps = grouped_peptides[pep_key];
        if (peps.length < 1) {
            return;
        }
        var first_pep = peps[0];
        var quant = null;
        var high_sn = false;
        var hexnac_type = {};
        peps.forEach(function(pep) {
            if ("CalculatedRatio" in pep) {
                quant = pep.CalculatedRatio;
                high_sn = true;
            }
            if ("has_low_sn" in pep && ! pep.has_low_sn) {
                high_sn = true;
            }
            if ("has_pair" in pep && pep.has_pair === true && pep.activation !== 'HCD') {
                quant = pep.QuanChannelID[0] == 1 ? 'potential_light' : 'potential_medium';
            }
            if ("hexnac_type" in pep) {
                hexnac_type[pep.hexnac_type] = true;
                if (("GlcNAc" in hexnac_type || "GalNAc" in hexnac_type)) {
                    delete hexnac_type['Unknown'];
                }
            }
        });
        var block = {
            'multi_protein' : false,
            'sequence' : first_pep.Sequence
        };
        if (first_pep.uniprot.length > 1) {
            block.multi_protein = true;
        }
        if (quant !== null && high_sn) {
            block.quant = quant;
        }
        if (Object.keys(hexnac_type).length > 0) {
            block.hexnac_type = Object.keys(hexnac_type);
        }
        if (first_pep.modifications) {
            block.sites = first_pep.modifications;
        }
        if (first_pep.Composition) {
            block.composition = first_pep.Composition;
        } else if (first_pep.modifications) {
            var count = 0;
            block.composition = [ [first_pep.modifications.map(function(site) { count += 1; return site[1]; })[0], count ].reverse().join('x') ];
        }
        first_pep.uniprot.forEach(function(uniprot) {
            if ( ! data[uniprot] ) {
                data[uniprot] = [];
            }
            data[uniprot].push(block);
        });
    });
    return data;
};


var partial = function(fn) {
    var args = Array.prototype.slice.call(arguments, 1);
    return function() { return fn.apply(this, args.concat(Array.prototype.slice.call(arguments, 0))); };
};

var onlyUnique = function(value, index, self) {
    return self.indexOf(value) === index;
};

var pbcopy = function(data) { var proc = require('child_process').spawn('pbcopy'); proc.stdin.write(data); proc.stdin.end(); }

var print_pep = function(protein,pep) { return [ protein, pep.sequence, pep.quant, pep.composition.map(function(comp) { return comp.replace(/x.*/,''); }).join(',') ].join('\t'); };

var filter_global_results = function(datablock,filter) {
if ( typeof datablock === 'function' ) {
    filter = datablock;
    datablock = global_datablock;
}
if ( ! filter ) {
    filter = function() { return true; };
}
return "uniprot\tsequence\tquant\tglyco\n"+
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
}


var global_results = [];
var global_datablock;

var db = new sqlite3.Database(files_to_open[0],sqlite3.OPEN_READONLY,function(err) {

    promisify_sqlite(db);
    init_caches(db).then(function() {
        console.log("Searching Peptides");
        global_db = db;

        // Area is the sum of Areas for the Events asociated with the QuanResultID (from EventAnnotations)
        retrieve_quantified_peptides(db).then(partial(check_quantified_peptides,db)).then(function(peps) {
            db.end_statement(peptide_metadata_sql);
            console.log("Done singlet checking");
            global_results = global_results.concat(peps);
        });
        Promise.reject(false).then(partial(retrieve_ambiguous_peptides,db)).then(function(peps) {
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
                if (nconf.get('hcd-processing-node')) {
                    processing_node = parseInt(nconf.get('hcd-processing-node'));
                }
                if (to_cut.length > 0) {
                    return Promise.all(  to_cut.splice(0,split_length).map(function(pep) { return get_spectrum(db,pep,processing_node).then( partial(check_galnac_glcnac_ratio,pep) ); }) ).then(arguments.callee);
                }
                return peps;
            });
            return result;
        }).then(function(peps) {
            console.log("Got all data");
            global_datablock = combine_all_peptides(peps);
            if (nconf.get('write-peptides')) {
                process.stdout.write(filter_global_results(function(pep) { return true; }));
                gui.App.quit();
            }
        }).catch(function(err) { console.log(arguments); console.log("Rejected"); });

    });


});