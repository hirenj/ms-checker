var combiner = require('./combiner');

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

var produce_peptide_data = function(db,pep) {
    return db.do_statement(peptide_metadata_sql, [ pep.PeptideID ]).then(function(pep_datas) {
        if ( ! pep.uniprot ) {
            pep.uniprot = [];
        }
        pep_datas.forEach(function(pep_data) {
            var uniprot = pep_data.Description.split('|')[1];
            if (uniprot && pep.uniprot.indexOf(uniprot) < 0) {
                pep.uniprot.push(uniprot);
            }
            if ( ! pep.pep_start ) {
                pep.pep_start = pep_data.Sequence.indexOf(pep.Sequence);
            }
        });
    });
};

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
        pep.modifications.forEach(function(mod) {
            mod[0] = pep.pep_start + mod[0];
        });
    }
    if (pep.possible_mods) {
        pep.possible_mods.forEach(function(mod) {
            mod[0] = pep.pep_start + mod[0];
            if (mod[3]) {
                mod[3] = pep.pep_start + mod[3];
                mod[4] = pep.pep_start + mod[4];
            }
        });
    }
};

exports.produce_peptide_data = produce_peptide_data;
exports.cleanup = function(db) {
    db.end_statement(peptide_metadata_sql);
};
exports.composition = extract_composition;
exports.combine = combiner.combine;
exports.modification_key = modification_key;
exports.fix_site_numbers = fix_ids;