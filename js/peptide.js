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
                peptide_scores_cache[score.PeptideID] = score.ScoreValue;
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
        if (peps.map(function(pep) { return pep.Sequence; }).filter(onlyUnique).length > 1) {
            bad_spectra.push(parseInt(spectrumID));
        }
    });
    return all_peps.filter(function(pep) { return bad_spectra.indexOf(pep.SpectrumID) < 0; });
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
exports.modification_key = modification_key;
exports.fix_site_numbers = fix_ids;
exports.filter_ambiguous_spectra = filter_ambiguous_spectra;