var util = require('util');
var peptide_search = require('./peptide.js');
var nconf = require('nconf');

var mod_string = peptide_search.mod_string;

Math.median = require('../js/math-median').median;

var Quantitative = function Quantitative() {

};

util.inherits(Quantitative,require('./processing-step.js'));

module.exports = exports = new Quantitative();


var wt_channel = "light";
var ko_channel = "medium";

exports.swapChannelLabels = function(do_swap) {
    if (do_swap) {
        wt_channel = "medium";
        ko_channel = "light";
    } else {
        wt_channel = "light";
        ko_channel = "medium";
    }
}

var dimethyl_names = {
    'Dimethyl:2H4' : 'medium',
    'Dimethyl' : 'light'
};

Object.defineProperty(exports, "wt_channel", { get: function () { return wt_channel; } });
Object.defineProperty(exports, "ko_channel", { get: function () { return ko_channel; } });

const dimethyl_counts_sql = 'SELECT \
    count(distinct Position) as count, PeptideID \
FROM PeptidesAminoAcidModifications \
    JOIN AminoAcidModifications USING(AminoAcidModificationID) \
WHERE (ModificationName = "Dimethyl" OR ModificationName = "Dimethyl:2H(4)") GROUP BY PeptideID';

const dimethyl_total_count_sql = 'SELECT count(distinct PeptideID) as count FROM PeptidesAminoAcidModifications';

const related_quants_sql = 'SELECT \
    EventID, \
    FileID, \
    Mass, \
    RT, \
    LeftRT, \
    RightRT, \
    SN, \
    Intensity, \
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


// We don't use this query because the MSF
// is missing the index on this table, and
// so is dog slow

const peptide_modification_sql = 'SELECT \
    Position,AminoAcidModificationID \
FROM PeptidesAminoAcidModifications \
WHERE PeptideID = ?';


const quant_channel_config_sql = 'SELECT \
    ParameterValue \
FROM ProcessingNodeParameters \
WHERE ParameterName = "QuantificationMethod"';


const glyco_peptide_modifications_sql = 'SELECT \
    PeptideID, Position, ModificationName, DeltaMass \
FROM PeptidesAminoAcidModifications \
    LEFT JOIN AminoAcidModifications USING (AminoAcidModificationID) \
WHERE \
    PeptideID in (SELECT distinct PeptideID \
    FROM PeptidesAminoAcidModifications \
    WHERE AminoAcidModificationID in (SELECT AminoAcidModificationID \
        FROM AminoAcidModifications \
        WHERE ModificationName like "%'+mod_string+'%") ) \
    OR PeptideID in (SELECT PeptideID \
        FROM \
            FileInfos LEFT JOIN MassPeaks USING(FileID) \
            LEFT JOIN SpectrumHeaders USING(MassPeakID) \
            JOIN Peptides using(SpectrumID) \
        WHERE FileInfos.Filename NOT LIKE "%raw" \
            AND Peptides.ConfidenceLevel = 3 \
            ) \
UNION \
\
SELECT \
    PeptideID, -1, ModificationName, DeltaMass \
FROM PeptidesTerminalModifications \
    LEFT JOIN AminoAcidModifications ON PeptidesTerminalModifications.TerminalModificationID = AminoAcidModifications.AminoAcidModificationID \
WHERE \
    PeptideID in (SELECT distinct PeptideID \
    FROM PeptidesAminoAcidModifications \
    WHERE AminoAcidModificationID in (SELECT AminoAcidModificationID \
        FROM AminoAcidModifications \
        WHERE ModificationName like "%'+mod_string+'%") ) \
    OR PeptideID in (SELECT PeptideID \
        FROM \
            FileInfos LEFT JOIN MassPeaks USING(FileID) \
            LEFT JOIN SpectrumHeaders USING(MassPeakID) \
            JOIN Peptides using(SpectrumID) \
        WHERE FileInfos.Filename NOT LIKE "%raw" \
            AND Peptides.ConfidenceLevel = 3 \
            ) \
';

const all_peptide_modifications_sql = 'SELECT \
    PeptideID, Position, ModificationName, DeltaMass \
FROM PeptidesAminoAcidModifications \
    LEFT JOIN AminoAcidModifications USING (AminoAcidModificationID) \
UNION \
\
SELECT \
    PeptideID, -1, ModificationName, DeltaMass \
FROM PeptidesTerminalModifications \
    LEFT JOIN AminoAcidModifications ON PeptidesTerminalModifications.TerminalModificationID = AminoAcidModifications.AminoAcidModificationID \
';

const all_peptide_modifications_count_sql = 'SELECT \
    count(*) as count \
FROM PeptidesAminoAcidModifications \
UNION \
SELECT \
    count(*) as count \
FROM PeptidesTerminalModifications \
';

const search_peptides_count_sql = 'SELECT count(*) as count from Peptides where ConfidenceLevel = 3 and SearchEngineRank <= 4';

// We need to merge with ReporterIonQuanResults when this table exists
// LEFT JOIN ReporterIonQuanResults USING(SpectrumID)

// SELECT count(*) FROM sqlite_master WHERE type='table' AND name='ReporterIonQuanResults'

const search_peptides_sql = 'SELECT \
    peptides.PeptideID, \
    peptides.UniquePeptideSequenceID, \
    peptides.SpectrumID, \
    peptides.Sequence, \
    PrecursorIonQuanResultsSearchSpectra.QuanResultID, \
    PrecursorIonQuanResults.QuanChannelID, \
    PrecursorIonQuanResults.Area, \
    ScanEvents.ActivationType as ActivationType, \
    peptides.SearchEngineRank, \
    SpectrumHeaders.Mass as mass, \
    SpectrumHeaders.Charge as charge, \
    MassPeaks.FileID as FileID, \
    ifnull(customdata.quantcount,0) as acceptedquant \
FROM peptides \
    LEFT JOIN PrecursorIonQuanResultsSearchSpectra \
        ON peptides.SpectrumID = PrecursorIonQuanResultsSearchSpectra.SearchSpectrumID \
    LEFT JOIN SpectrumHeaders USING(SpectrumID) \
    LEFT JOIN ScanEvents USING(ScanEventID) \
    LEFT JOIN MassPeaks USING(MassPeakID) \
    LEFT JOIN PrecursorIonQuanResults \
        ON PrecursorIonQuanResultsSearchSpectra.QuanResultID = PrecursorIonQuanResults.QuanResultID \
    LEFT JOIN (SELECT PeptideID, count(FieldID) as quantcount FROM CustomDataPeptides \
                WHERE FieldID in (SELECT FieldID FROM CustomDataFields WHERE DisplayName = "QuanResultID") GROUP BY PeptideID) as customdata \
        ON peptides.PeptideID = customdata.PeptideID  \
WHERE peptides.ConfidenceLevel = 3 \
AND SearchEngineRank <= 4 \
AND peptides.PeptideID in (SELECT distinct PeptideID \
    FROM PeptidesAminoAcidModifications \
    WHERE AminoAcidModificationID in (SELECT AminoAcidModificationID \
        FROM AminoAcidModifications \
        WHERE ModificationName like "%'+mod_string+'%") \
    )';

const search_peptides_sql_tmt = 'SELECT \
    peptides.PeptideID, \
    peptides.UniquePeptideSequenceID, \
    peptides.SpectrumID, \
    peptides.Sequence, \
    ReporterIonQuanResults.QuanChannelID, \
    ReporterIonQuanResults.Height, \
    ScanEvents.ActivationType as ActivationType, \
    peptides.SearchEngineRank, \
    SpectrumHeaders.Mass as mass, \
    SpectrumHeaders.Charge as charge, \
    MassPeaks.FileID as FileID \
FROM peptides \
    LEFT JOIN SpectrumHeaders USING(SpectrumID) \
    LEFT JOIN ScanEvents USING(ScanEventID) \
    LEFT JOIN MassPeaks USING(MassPeakID) \
    LEFT JOIN ReporterIonQuanResultsSearchSpectra \
    ON peptides.SpectrumID = ReporterIonQuanResultsSearchSpectra.SearchSpectrumID \
    LEFT JOIN ReporterIonQuanResults \
    ON ReporterIonQuanResultsSearchSpectra.SpectrumID = ReporterIonQuanResults.SpectrumID \
WHERE peptides.ConfidenceLevel = 3 \
AND SearchEngineRank <= 4 \
AND peptides.PeptideID in (SELECT distinct PeptideID \
    FROM PeptidesAminoAcidModifications \
    WHERE AminoAcidModificationID in (SELECT AminoAcidModificationID \
        FROM AminoAcidModifications \
        WHERE ModificationName like "%'+mod_string+'%") \
    )';

const MASS_MEDIUM = 32.056407;
const MASS_LIGHT = 28.0313;

var dimethyl_count_cache = null;

var find_dimethyls = function(db,pep) {
    return new Promise(function(resolve,reject) {

        if (! dimethyl_count_cache ) {

            dimethyl_count_cache = {};
            exports.notify_task('Populating Dimethyl cache');
            exports.notify_progress(0,1);
            dimethyl_count_cache = {};
            var total_count = 1000000;
            var idx = 0;
            db.all(dimethyl_total_count_sql).then(function(count) {
                total_count = count[0].count;
            });
            db.each(dimethyl_counts_sql,[],function(err,count) {
                idx += 1;
                exports.notify_progress(idx,total_count);
                dimethyl_count_cache[count.PeptideID] = count.count+1;
            }).then(function(data) {
                exports.notify_progress(total_count,total_count);

                if ( ! pep.QuanResultID ) {
                    resolve();
                } else {
                    resolve(find_dimethyls(db,pep));
                }
            },reject);

        } else {
            resolve(dimethyl_count_cache[pep.PeptideID]);
        }
    });
};


var validated_quans_cache = {};


/*
FIXME: Write Unit tests for this section - make up some minimal MSF with the peaks that
we need to find, and one where there is a singlet. Possible minimal method call where
we can just grab the pair for a given peptide too?
*/

var check_potential_pair = function(db,pep,num_dimethyl,channel_conf) {
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
                resolve(check_potential_pair(db,pep,count,channel_conf));
            },reject);
            return;
        }

        pep.has_pair = false;

        var mass_change_dir = 1;
        if (pep.QuanChannelID.indexOf("medium") >= 0) {
            mass_change_dir = -1;
        }
        var sns = [];
        var intensities = [];
        pep.has_low_sn = false;
        pep.has_high_sn = false;
        var numeric_channel_id = channel_conf[pep.QuanChannelID.filter(onlyUnique)[0]];

        db.all(related_quants_sql,[ pep.QuanResultID, numeric_channel_id ]).then(function(events) {
            var events_length = events.length;
            events.forEach(function(ev) {

                sns.push(ev.SN);
                intensities.push(ev.Intensity);

                events_length -= 1;
                if (events_length <= 0) {

                    pep.quant_intensity_range = [ Math.min.apply(null,intensities), Math.max.apply(null,intensities) ];

                    if (Math.median(sns) >= 100 ) {
                        pep.has_high_sn = true;
                    }
                    if (Math.median(sns) < 10) {
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
                        // if (nconf.get('MS_DEBUG') || nconf.get('debug')) {
                        //     pep.matching_events = rows.map(function(row) { return row.EventID; });
                        //     pep.search_event = ev.EventID;
                        //     pep.target_mass = target_mass;
                        //     pep.num_dimethyl = num_dimethyl;
                        //     pep.mass_change_dir = mass_change_dir;
                        // }
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


// Channel identifier found in
// ProcessingNodeParameters.ParameterName QuantificationMethod

// We wish to combine peptides that have the same Peptide ID, but have
// different quant channels so that we don't have repeated rows for different
// components of the same quantification result.

var combine_quantified_peptides = function(peps,channel_conf) {
    var all_peps = {};
    peps.forEach(function(pep) {
        var curr_pep = all_peps[pep.PeptideID] ||  pep;
        curr_pep.areas = curr_pep.areas || {};
        curr_pep.heights = curr_pep.heights || {};

        if ( pep.QuanChannelID && pep.QuanChannelID != "medium" || pep.QuanChannelID != "light" ) {
            pep.QuanChannelID = channel_conf[pep.QuanChannelID];
        }

        if (pep.Area && pep.QuanChannelID) {
            if (pep.QuanChannelID != 'medium' && pep.QuanChannelID != 'light') {
                throw new Error("Could not resolve Quant channel ID",pep.QuanChannelID);
            }
            curr_pep.areas[pep.QuanChannelID] = pep.Area;
        }
        if (pep.Height != null && pep.QuanChannelID) {
            curr_pep.heights[pep.QuanChannelID] = pep.Height;
        }
        if (! Array.isArray(curr_pep.QuanChannelID) ) {
            curr_pep.QuanChannelID =  curr_pep.QuanChannelID ? [ curr_pep.QuanChannelID ] : [];
        }
        if (curr_pep != pep && pep.QuanChannelID) {
            curr_pep.QuanChannelID.push(pep.QuanChannelID);
        }
        if (Object.keys(curr_pep.areas).length != curr_pep.QuanChannelID.length && Object.keys(curr_pep.heights).length != curr_pep.QuanChannelID.length) {
            throw new Error("quant_channel_counts")
        }
        all_peps[pep.PeptideID] = curr_pep;
    });
    return Object.keys(all_peps).map(function(key) { return all_peps[key]; });
};

var extract_channel_info = function(channel_info) {
    var params = channel_info.Parameter.map(function(param) {
        var obj = {};
        obj[ param['$'].name ] =  param['_'];
        return obj;
    });
    params = params.reduce(function(prev,curr) {
        for (var key in curr) {
            prev[key] = curr[key];
        }
        return prev;
    });
    return params;
};

var get_channel_config = function(db) {
    return db.all(quant_channel_config_sql).then(function(xml) {
        return new Promise(function(resolve,reject) {
            if ( ! xml || xml.length < 1 ) {
                resolve([]);
            }
            require('xml2js').parseString(xml[0].ParameterValue,function(err,result) {
                if (err) {
                    throw err;
                }
                var channel_infos = result.ProcessingMethod.MethodPart.filter(function(part) {
                    return part.$.name == 'QuanChannels' || part.$.name == 'MassTags';
                })[0].MethodPart;
                var extracted = channel_infos.map(extract_channel_info);
                resolve(extracted);
            });
        });
    });
}

var retrieve_quantified_peptides = function(db) {
    var peps = [];
    exports.notify_task('Retrieving quantified peptides');
    exports.notify_progress(0,1);
    var total_count = 100000;
    var idx = 0;
    db.all(search_peptides_count_sql).then(function(count) {
        total_count = count[0].count;
    });

    let is_tmt = true;

    let search_sql = is_tmt ? search_peptides_sql_tmt : search_peptides_sql;

    return db.each(search_sql,[],function(err,pep) {
        idx += 1;
        if (pep) {
            pep.acceptedquant = (pep.acceptedquant == 0) ? false : true;
            peps.push(pep);
        }
        exports.notify_progress(idx,total_count);
    }).then(function() {
        exports.notify_progress(total_count,total_count);
        return peps;
    });
};

var onlyUnique = function(value, index, self) {
    return self.indexOf(value) === index;
};

var check_quantified_peptides = function(db,peps) {
    var idx = 0;
    var update_fractions = function(pep) {
        idx += 1;
        exports.notify_progress(idx,total_peps);
    };
    var setup_config = get_channel_config(db).then(function(quant_config) {
        var channel_conf = {};
        if (quant_config.length != 2) {
            // This only works because Dimethyl has QuanResultID
            var quant_peps = peps.filter(function(pep) { return pep.QuanResultID; });
            if (quant_peps.length > 0) {
                throw new Error("Incorrect number of quant channels");
            }
        }
        quant_config.forEach(function(conf) {
            if (conf.ChannelName) {
                channel_conf[ conf.ChannelName.toLowerCase() ] = conf.ChannelID;
                channel_conf[ conf.ChannelID ] = conf.ChannelName.toLowerCase();
            }
            if (conf.TagName) {
                channel_conf[ conf.TagName.toLowerCase() ] = conf.TagID;
                channel_conf[ conf.TagID ] = conf.TagName.toLowerCase();                
            }
        });
        return channel_conf;
    });
    var combined_peps = [];
    return setup_config.then(function(channel_conf) {
        combined_peps = combine_quantified_peptides(peps,channel_conf);

        var singlet_peps = combined_peps.filter( function(pep) { return pep.QuanChannelID.filter(onlyUnique).length == 1; } );
        var total_peps = singlet_peps.length + combined_peps.length;

        exports.notify_task('Validating quantified peptides');
        exports.notify_progress(0,total_peps);

        var pair_promises = singlet_peps.map( function(pep) { return check_potential_pair(db,pep,null,channel_conf).then(update_fractions); } );
        var metadata_promises = combined_peps.map( function(pep) { return peptide_search.produce_peptide_data(db,pep).then( produce_peptide_modification_data(db,pep) ).then(update_fractions);  } );

        return pair_promises.concat(metadata_promises);
    }).then(function(promises) {
        return Promise.all(promises);
    }).then(function() {
        peptide_search.cleanup(db);
        exports.notify_progress(total_peps,total_peps);
        return combined_peps;
    });
};

var peptide_modifications_cache = {'empty' : true };

var clean_mod = function(mod) {
    var new_mod = [].concat(mod);
    new_mod[1] = new_mod[1].replace(/\d/g,'');
    return new_mod;
};

var produce_peptide_modification_data = function(db,pep) {
    return new Promise(function(resolve,reject) {

        if (! peptide_modifications_cache || peptide_modifications_cache.empty ) {

            delete peptide_modifications_cache['empty'];

            exports.notify_task('Population of modifications cache');
            exports.notify_progress(0,1);

            var total_count = 100000000;
            var idx = 0;

            db.all(all_peptide_modifications_count_sql).then(function(rows) {
                if (rows.length > 1) {
                    total_count = rows[0].count + rows[1].count;
                }
            });
            db.each(glyco_peptide_modifications_sql,[],function(err,mods) {
                idx += 1;
                peptide_modifications_cache[mods.PeptideID] =  peptide_modifications_cache[mods.PeptideID] || [];
                var mod_name = mods.ModificationName.replace(/-/g,'');
                mod_name = mod_name.replace(/[\(\)]/g,'');
                peptide_modifications_cache[mods.PeptideID].push([ mods.Position < 0 ? 1 : mods.Position + 1, mod_name, mods.DeltaMass ]);
                exports.notify_progress(idx,total_count);
            }).then(function() {
                exports.notify_progress(total_count,total_count);
                if ( ! pep.PeptideID ) {
                    resolve();
                } else {
                    resolve(produce_peptide_modification_data(db,pep));
                }
            },reject);
        } else {
            var dimethyls = (peptide_modifications_cache[pep.PeptideID] || []).map( mod => dimethyl_names[mod[1]] ).filter( mod => mod ).filter(onlyUnique);
            pep.dimethyl_modification = dimethyls[0];
            pep.modifications = (peptide_modifications_cache[pep.PeptideID] || []).filter( function(mod) { return (mod[1] || '').indexOf(mod_string) >= 0 } ).map( clean_mod );
            resolve();
        }
    });
};

var populate_missing_quant_channel = function(peps) {
    peps.forEach( pep => {
        if (! pep.dimethyl_modification) {
            var dimethyls = (peptide_modifications_cache[pep.PeptideID] || []).map( mod => dimethyl_names[mod[1]] ).filter( mod => mod ).filter(onlyUnique);
            pep.dimethyl_modification = dimethyls[0];
        }
    });
    return peps;
};

var init_caches = function(db) {
    return produce_peptide_modification_data(db,{}).then(function() { return find_dimethyls(db,{'PeptideID': 0}); });
};

var clear_caches = function() {
    validated_quans_cache = {};
    dimethyl_count_cache = null;
    peptide_modifications_cache = {'empty' : true };
    exports.modifications_cache = peptide_modifications_cache;
};

exports.init_caches = init_caches;
exports.clear_caches = clear_caches;
exports.retrieve_quantified_peptides = retrieve_quantified_peptides;
exports.modifications_cache = peptide_modifications_cache;
exports.populate_missing_quant_channel = populate_missing_quant_channel;
exports.check_quantified_peptides = check_quantified_peptides;
