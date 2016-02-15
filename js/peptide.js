var util = require('util');

const peptide_metadata_sql = 'SELECT \
    Description, \
    Proteins.Sequence \
FROM Peptides \
    LEFT JOIN PeptidesProteins USING(PeptideID) \
    LEFT JOIN ProteinAnnotations USING (ProteinID) \
    LEFT JOIN Proteins USING (ProteinID) \
WHERE Peptides.PeptideID = ?';

const protein_sequence_sql = 'SELECT \
    Sequence \
FROM Proteins \
WHERE ProteinID = ?'

const peptide_score_sql = 'SELECT \
    PeptideID, \
    ScoreValue \
FROM PeptideScores';

const peptide_score_count_sql = 'SELECT count(*) as count from PeptideScores';

var Peptide = function Peptide() {

};

util.inherits(Peptide,require('./processing-step.js'));

module.exports = exports = new Peptide();

var onlyUnique = function(value, index, self) {
    return self.indexOf(value) === index;
};

var peptide_scores_cache = {'empty' : true};

var produce_peptide_score_data = function(db,pep) {
    var self = this;
    return new Promise(function(resolve,reject) {

        if (! peptide_scores_cache || peptide_scores_cache.empty ) {

            delete peptide_scores_cache['empty'];

            exports.notify_task('Population of peptide scores');
            exports.notify_progress(0,1);

            var total_count = 100000000;
            var idx = 0;

            db.all(peptide_score_count_sql).then(function(rows) {
                if (rows.length > 1) {
                    total_count = rows[0].count;
                }
            });
            db.each(peptide_score_sql,[],function(err,score) {
                if ( score ) {
                    peptide_scores_cache[score.PeptideID] = score.ScoreValue;
                }
                idx += 1;
                exports.notify_progress(idx,total_count);
            }).then(function() {
                exports.notify_progress(total_count,total_count);
                if ( ! pep.PeptideID ) {
                    resolve();
                } else {
                    resolve(produce_peptide_score_data(db,pep));
                }
            });
        } else {
            resolve(peptide_scores_cache[pep.PeptideID]);
        }
    });
};


var produce_peptide_data = function(db,pep) {
    return db.do_statement(peptide_metadata_sql, [ pep.PeptideID ]).then(function(pep_datas) {
        if ( ! pep.uniprot ) {
            pep.uniprot = [];
        }
        if ( ! pep.gene ) {
            pep.gene = [];
        }
        var starts = [];
        pep_datas.forEach(function(pep_data) {
            if ( ! pep_data.Description ) {
                return;
            }
            var uniprot = pep_data.Description.split('|')[1];
            if (uniprot && pep.uniprot.indexOf(uniprot) < 0) {
                pep.uniprot.push(uniprot);
            }
            var genes;
            var re = /(?:GN=)([^\s]+)/g;
            while( genes = re.exec(pep_data.Description) ) {
                pep.gene.push(genes[1]);
            }
            starts.push(pep_data.Sequence.indexOf(pep.Sequence));

            if ( ! pep.pep_start ) {
                pep.pep_start = pep_data.Sequence.indexOf(pep.Sequence);
            }
        });
        if (starts.length > 1) {
            pep.starts = starts;
        }
    });
};

var produce_peptide_scores = function(db,peps) {
    return Promise.all(peps.map(function(pep) {
        return produce_peptide_score_data(db,pep).then(function(score) {
            pep.score = score;
            return pep;
        });
    }));
};

var filter_ambiguous_spectra = function(all_peps) {
    var peptides_by_spectrum = {};
    var bad_spectra = [];
    all_peps.forEach(function(pep) {
        if ( ! peptides_by_spectrum[pep.SpectrumID]) {
            peptides_by_spectrum[pep.SpectrumID] = [];
        }
        peptides_by_spectrum[pep.SpectrumID].push(pep);
    });
    Object.keys(peptides_by_spectrum).forEach(function(spectrumID) {
        var peps = peptides_by_spectrum[spectrumID];
        var seqs = peps.map(function(pep) { return pep.Sequence; }).filter(onlyUnique);
        var lengths = seqs.map(function(seq) { return seq.length; }).filter(onlyUnique);

        // If there's only one sequence, we don't have to do anything

        if (seqs.length < 2) {
            return;
        }

        // If we have more than one length of sequence mapped to this spectrum
        // we don't want to keep it

        if (lengths.length > 1) {
            bad_spectra.push(parseInt(spectrumID));
            return;
        }

        // We need to look for single amino acid differences between peptides
        // to look for the Leucine / Isoleucine difference between peptides.

        var aas = seqs.map(function(seq) {  return seq.split(''); });
        var diff_string = '';
        for (var i = 0; i < aas[0].length; i++) {
            var pos_aas = aas.map(function(aa_set) { return aa_set[i]; }).filter(onlyUnique);
            if (pos_aas.length < 2) {
                return;
            }
            diff_string = diff_string + pos_aas.join('');
        }

        // If the difference between the two sequences only consists of I and L
        // then the sequences are the same.

        if (diff_string.match(/^[IL]+$/)) {
            return;
        }


        // If we don't match the Leucine / Isoleucine get out of jail
        // free card, we should mark this as a bad spectrum

        bad_spectra.push(parseInt(spectrumID));
    });
    return all_peps.filter(function(pep) { return bad_spectra.indexOf(pep.SpectrumID) < 0; });
};

var calculate_deltacn = function(all_peps) {
    var grouped_peptides = {};

    all_peps.forEach(function(pep) {
        var key = pep.SpectrumID;
        if ( ! grouped_peptides[key] )  {
            grouped_peptides[key] = [];
        }
        grouped_peptides[key].push(pep);
    });
    Object.keys(grouped_peptides).forEach(function(pep_id) {
        var peps = grouped_peptides[pep_id].sort(function(a,b) {
            return a.score - b.score;
        }).reverse();
        if (peps.length > 1) {
            var last_score = peps[0].score;
            peps.forEach(function(pep) {
                pep.deltacn = (last_score - pep.score) / last_score;
            });
        } else {
            peps[0].deltacn = 0;
        }
    });
    return all_peps;
};

var filter_deltacn = function(cutoff,all_peps) {
    return all_peps.filter(function(pep) { return pep.deltacn <= cutoff; });
};

var merge_modifications_deltacn = function(all_peps) {
    var grouped_peptides = {};

    all_peps.forEach(function(pep) {
        var key = pep.SpectrumID;
        if ( ! grouped_peptides[key] )  {
            grouped_peptides[key] = [];
        }
        grouped_peptides[key].push(pep);
    });
    Object.keys(grouped_peptides).forEach(function(pep_id) {
        var peps = grouped_peptides[pep_id];
        if (peps.length > 1) {
            // We should split the peptides into ambiguous and unambiguous
            var ambiguous = peps.filter(function(pep) { return pep.possible_mods; });
            var unambiguous = peps.filter(function(pep) { return pep.modifications; });

            // If there is only a single unambiguous, then we just use the unambiguous and drop the other peptides
            if (unambiguous.length < 2) {
                ambiguous.forEach(function(pep) { pep.drop = true; });
                return;
            }

            // More than one unambiguous should merge together to form an ambiguous
            if (unambiguous.length > 1) {
                ambiguous.push(merge_unambiguous(unambiguous));
            }

            // All the ambiguous should then be merged together, placed onto the first peptide
            // and then the other peptides should be dropped.
            merge_ambiguous(ambiguous);
            debugger;
        }
    });
    return all_peps.filter(function(pep) { return ! pep.drop; });
};

var merge_unambiguous = function(peps) {
    debugger;
    var pep_positions = peps.map(function(pep,i) { if ( i > 0 ) { pep.drop = true; } return pep.modifications[0][0]; });

    // This should really be merged per composition
    var min_pos = Math.min.apply(Math,pep_positions);
    var max_pos = Math.max.apply(Math,pep_positions);
    peps[0].modifications.forEach(function(mod) {
        mod[3] = min_pos;
        mod[4] = max_pos;
    });
    peps[0].made_ambiguous = true;
    peps[0].possible_mods = peps[0].modifications;
    debugger;
    return peps[0];
}

var merge_ambiguous = function(peps) {
    return;
}

var extract_composition = function(mods) {
    var compositions = {};
    mods.forEach(function(mod) {
        if (! compositions[mod[1]]) {
            compositions[mod[1]] = 0;
        }
        compositions[mod[1]] += 1;
    });
    return Object.keys(compositions).map(function(comp) {  return compositions[comp]+"x"+comp; }).sort();
};

var modification_key = function(pep) {
    var mod_key = null;
    if (pep.modifications) {
        mod_key = pep.modifications.sort(function(a,b) { return a[0] - b[0]; }).map(function(mod) {  return mod[0]+"-"+mod[1]; }).join(',');
    }
    if (pep.Composition) {
        mod_key = pep.Composition.sort().join(',');
    }
    return mod_key;
};

var fix_ids = function(pep) {
    if (pep.modifications) {
        pep.modifications = pep.modifications.map(function(mod) { return [].concat(mod); });
        pep.modifications.forEach(function(mod) {
            mod[0] = pep.pep_start + mod[0];
        });
    }
    if (pep.possible_mods) {
        pep.possible_mods = pep.possible_mods.map(function(mod) { return [].concat(mod); });
        pep.possible_mods.forEach(function(mod) {
            mod[0] = pep.pep_start + mod[0];
            if (mod[3]) {
                mod[3] = pep.pep_start + mod[3];
                mod[4] = pep.pep_start + mod[4];
            }
        });
    }
};

var init_caches = function(db) {
    return produce_peptide_score_data(db,{});
};

exports.init_caches = init_caches;
exports.produce_peptide_data = produce_peptide_data;
exports.produce_peptide_scores = produce_peptide_scores;
exports.produce_peptide_scores_and_cleanup = function(db,peps) {
    return init_caches(db).then(function() {
        return produce_peptide_scores(db,peps);
    }).then(function(peps) {
        peptide_scores_cache = {'empty' : true};
        return peps;
    });
};
exports.cleanup = function(db) {
    db.end_statement(peptide_metadata_sql);
};
exports.composition = extract_composition;
exports.combine = function(peps) {
    return require('./combiner').combine(peps);
}
exports.calculate_deltacn = calculate_deltacn;
exports.filter_deltacn = filter_deltacn;
exports.merge_modifications_deltacn = merge_modifications_deltacn;
exports.modification_key = modification_key;
exports.fix_site_numbers = fix_ids;
exports.filter_ambiguous_spectra = filter_ambiguous_spectra;