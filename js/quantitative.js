var util = require('util');
var peptide_search = require('./peptide.js');

Math.median = require('../js/math-median').median;

var Quantitative = function Quantitative() {

};

util.inherits(Quantitative,require('events').EventEmitter);

module.exports = exports = new Quantitative();

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


const all_peptide_modifications_sql = 'SELECT \
    PeptideID, Position, ModificationName, DeltaMass \
FROM PeptidesAminoAcidModifications \
    LEFT JOIN AminoAcidModifications USING (AminoAcidModificationID) \
WHERE \
    PeptideID in (SELECT distinct PeptideID \
    FROM PeptidesAminoAcidModifications \
    WHERE AminoAcidModificationID in (SELECT AminoAcidModificationID \
        FROM AminoAcidModifications \
        WHERE ModificationName like "%Hex%") ) \
UNION \
\
SELECT \
    PeptideID, 0, ModificationName, DeltaMass \
FROM PeptidesTerminalModifications \
    LEFT JOIN AminoAcidModifications ON PeptidesTerminalModifications.TerminalModificationID = AminoAcidModifications.AminoAcidModificationID \
WHERE \
    PeptideID in (SELECT distinct PeptideID \
    FROM PeptidesAminoAcidModifications \
    WHERE AminoAcidModificationID in (SELECT AminoAcidModificationID \
        FROM AminoAcidModifications \
        WHERE ModificationName like "%Hex%") ) \
';

const all_peptide_modifications_count_sql = 'SELECT \
    count(*) as count \
FROM PeptidesAminoAcidModifications \
UNION \
SELECT \
    count(*) as count \
FROM PeptidesTerminalModifications \
';

const search_peptides_count_sql = 'SELECT count(*) as count from Peptides where ConfidenceLevel = 3 and SearchEngineRank = 1';

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


const MASS_MEDIUM = 32.056407;
const MASS_LIGHT = 28.0313;

var dimethyl_count_cache = null;

var find_dimethyls = function(db,pep) {
    return new Promise(function(resolve,reject) {

        if (! dimethyl_count_cache ) {

            dimethyl_count_cache = {};
            exports.emit('task','Populating Dimethyl cache');
            exports.emit('progress',0);
            dimethyl_count_cache = {};
            var total_count = 1000000;
            var idx = 0;
            db.all(dimethyl_total_count_sql).then(function(count) {
                total_count = count[0].count;
            });
            var last_frac = 0;
            db.each(dimethyl_counts_sql,[],function(err,count) {
                idx += 1;
                frac = parseFloat( (idx / total_count).toFixed(2) );
                if (frac !== last_frac) {
                    last_frac = frac;
                    exports.emit('progress',frac );
                }
                dimethyl_count_cache[count.PeptideID] = count.count+1;
            }).then(function(data) {
                exports.emit('progress',1);

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
    var peps = [];
    exports.emit('task','Retrieving quantified peptides');
    exports.emit('progress',0);
    var last_frac = 0;
    var total_count = 100000;
    var idx = 0;
    db.all(search_peptides_count_sql).then(function(count) {
        total_count = count[0].count;
    });
    return db.each(search_peptides_sql,[],function(err,pep) {
        idx += 1;
        var frac = parseFloat((idx/total_count).toFixed(2));
        peps.push(pep);
        if (frac !== last_frac) {
            exports.emit('progress',frac);
            last_frac = frac;
        }
    }).then(function() {
        exports.emit('progress',1);
        return peps;
    });
};

var onlyUnique = function(value, index, self) {
    return self.indexOf(value) === index;
};

var check_quantified_peptides = function(db,peps) {
    var last_frac = 0;
    var idx = 0;
    var update_fractions = function(pep) {
        idx += 1;
        var frac = parseFloat((idx / total_peps).toFixed(2));
        if ( frac !== last_frac ) {
            exports.emit('progress',frac);
            last_frac = frac;
        }
    };

    var combined_peps = combine_quantified_peptides(peps);

    var singlet_peps = combined_peps.filter( function(pep) { return pep.QuanChannelID.filter(onlyUnique).length == 1; } );
    var total_peps = singlet_peps.length + combined_peps.length;

    exports.emit('task','Validating quantified peptides');
    exports.emit('progress',0);

    var pair_promises = singlet_peps.map( function(pep) { return check_potential_pair(db,pep,null).then(update_fractions); } );
    var metadata_promises = combined_peps.map( function(pep) { return peptide_search.produce_peptide_data(db,pep).then( produce_peptide_modification_data(db,pep) ).then(update_fractions);  } );

    return Promise.all(pair_promises.concat(metadata_promises)).then(function() {
        peptide_search.cleanup(db);
        exports.emit('progress',1);
        return combined_peps;
    }).catch(function(err) {
        console.log(err);
    });
};

var peptide_modifications_cache = {'empty' : true };

var produce_peptide_modification_data = function(db,pep) {
    return new Promise(function(resolve,reject) {

        if (! peptide_modifications_cache || peptide_modifications_cache.empty ) {

            peptide_modifications_cache = peptide_modifications_cache || {};

            delete peptide_modifications_cache['empty'];

            exports.emit('task','Population of modifications cache');
            exports.emit('progress',0);

            peptide_modifications_cache = {};
            var last_frac = 0;
            var total_count = 100000000;
            var idx = 0;

            db.all(all_peptide_modifications_count_sql).then(function(rows) {
                if (rows.length > 1) {
                    total_count = rows[0].count + rows[1].count;
                }
            });
            db.each(all_peptide_modifications_sql,[],function(err,mods) {
                idx += 1;
                var frac = parseFloat((idx/total_count).toFixed(2));
                peptide_modifications_cache[mods.PeptideID] =  peptide_modifications_cache[mods.PeptideID] || [];
                peptide_modifications_cache[mods.PeptideID].push([ mods.Position + 1, mods.ModificationName, mods.DeltaMass ]);
                if (frac !== last_frac) {
                    exports.emit('progress',frac);
                    last_frac = frac;
                }
            }).then(function() {
                exports.emit('progress',1);
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

exports.init_caches = init_caches;
exports.retrieve_quantified_peptides = retrieve_quantified_peptides;
exports.modifications_cache = peptide_modifications_cache;
exports.check_quantified_peptides = check_quantified_peptides;
