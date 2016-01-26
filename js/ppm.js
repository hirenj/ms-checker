var masses = require('./fragmentions').calculate_peptide_mass;
var quantitative = require('./quantitative');


// const retrieve_ppm_all_peptides_sql = 'SELECT \
//     Peptides.PeptideID as PeptideID, \
//     Peptides.Sequence as Sequence, \
//     SpectrumHeaders.Mass as Mass, \
//     PeptideScores.ScoreValue as ScoreValue, \
//     Position, \
//     ModificationName, \
//     DeltaMass \
// FROM Peptides LEFT JOIN  SpectrumHeaders ON (Peptides.SpectrumID = SpectrumHeaders.SpectrumID) LEFT JOIN PeptideScores ON (Peptides.PeptideID = PeptideScores.PeptideID) LEFT JOIN ( \
// SELECT \
//     PeptideID, Position, ModificationName, DeltaMass \
// FROM PeptidesAminoAcidModifications \
//     LEFT JOIN AminoAcidModifications USING (AminoAcidModificationID) \
// UNION \
// SELECT \
//     PeptideID, -1, ModificationName, DeltaMass \
// FROM PeptidesTerminalModifications \
//     LEFT JOIN AminoAcidModifications ON PeptidesTerminalModifications.TerminalModificationID = AminoAcidModifications.AminoAcidModificationID \
// ) as mods ON (Peptides.PeptideID = mods.PeptideID) \
// ';

const retrieve_ppm_all_peptides_sql = 'SELECT \
	Peptides.PeptideID as PeptideID, \
	Peptides.Sequence as Sequence, \
	SpectrumHeaders.Mass as Mass, \
	PeptideScores.ScoreValue as ScoreValue, \
	MassPeaks.FileID as FileID \
FROM Peptides LEFT JOIN  SpectrumHeaders USING(SpectrumID) LEFT JOIN MassPeaks USING(MassPeakID) LEFT JOIN PeptideScores USING(PeptideID) \
WHERE SpectrumHeaders.Charge BETWEEN 2 AND 4 \
';

const retrieve_mods_sql = 'SELECT DISTINCT \
	PeptideID, Position, ModificationName, DeltaMass \
FROM PeptidesAminoAcidModifications \
	LEFT JOIN AminoAcidModifications USING (AminoAcidModificationID) \
UNION \
SELECT \
	PeptideID, -1 as Position, ModificationName, DeltaMass \
FROM PeptidesTerminalModifications \
	LEFT JOIN AminoAcidModifications ON PeptidesTerminalModifications.TerminalModificationID = AminoAcidModifications.AminoAcidModificationID \
ORDER BY PeptideID';

// const retrieve_ppm_all_peptides_sql = 'SELECT \
//     Peptides.PeptideID as PeptideID, \
//     Peptides.Sequence as Sequence, \
//     SpectrumHeaders.Mass as Mass, \
//     PeptideScores.ScoreValue as ScoreValue \
// FROM Peptides LEFT JOIN  SpectrumHeaders ON (Peptides.SpectrumID = SpectrumHeaders.SpectrumID) LEFT JOIN PeptideScores ON (Peptides.PeptideID = PeptideScores.PeptideID) \
// ';

// We love Protein Prospector when trying to
// figure out what the right masses should be
// http://prospector.ucsf.edu/prospector/cgi-bin/mssearch.cgi

var calculate_theoretical_mass = function(peptide,modifications) {
	var aa_masses = masses(peptide,modifications);
	// We want to take the Mass column from SpectrumHeaders (which is the MH+) in the spreadsheets
	peptide.mass_theoretical = (aa_masses.reduce(function(p,n) { return p+n; },0));
};

var get_ppm = function(peptide,modifications) {
	calculate_theoretical_mass(peptide,modifications);
	peptide.ppm  = 1e06*((peptide.mass / peptide.mass_theoretical) - 1);
};

var annotate_ppm = function(peps) {
	peps.forEach(get_ppm);
	// foo = peps.map(function(pep) { return pep.ppm+"\t"+pep.score; }).join("\n");
	// require('fs').writeFileSync('foo.txt',foo);
	// process.exit(1);
	// debugger;
	return peps;
};

var derive_ppm_batch = function(pep) {
	get_ppm(pep,pep.modifications);
	return {'ppm' : pep.ppm, 'score' : pep.ScoreValue };
};

var retrieve_all_ppms = function(db) {
	var current_pep = { 'modifications' : [] };
	var all_mods = { 0 : [] };
	var all_values = {};

	var last_peptide_id = 0;

	var mods_function = function(err,mods) {
		if ( ! all_mods[mods.PeptideID] ) {
			all_mods[mods.PeptideID] = [];
		}
		var mod_name = mods.ModificationName.replace(/-/g,'');
		var position = mods.Position < 0 ? 1 : mods.Position + 1;
		// if (all_mods[mods.PeptideID].filter(function(mod) { return mod[0] == position && mod[1] == mod_name;  }).length == 0) {
			all_mods[mods.PeptideID].push([ position, mod_name, mods.DeltaMass ]);        	
		// }
		console.log(3246115 - mods.PeptideID);
		if (mods.PeptideID !== last_peptide_id) {
			all_mods[last_peptide_id] = all_mods[last_peptide_id].reduce(function(p,n) { return p + n[2] },0);
			last_peptide_id = mods.PeptideID;
		}
	};

	var all_mods_promise = db.each(retrieve_mods_sql,[],mods_function);

	//blah = by(foo$V3,cut(foo$V2,breaks=seq(-15,15,by=0.5)),function(x) {length(x[x > 3.104]) })
	//blah = by(foo$V3,cut(foo$V2,breaks=seq(-15,15,by=0.5)),function(x) { browser(); length(x[x > quantile(x,0.99)]) })

	return all_mods_promise.then(function() {
		return db.each(retrieve_ppm_all_peptides_sql,[],function(err,pep) {
			if ( ! all_values[pep.FileID]) {
				all_values[pep.FileID] = [];
			}
			console.log(3246115 - pep.PeptideID);
			if (pep.PeptideID !== current_pep.PeptideID && current_pep.PeptideID) {
				current_pep.modifications = [ [ 1, 1 , all_mods[current_pep.PeptideID] ] ];
				all_values[current_pep.FileID].push(derive_ppm_batch(current_pep));
				if (Math.abs(current_pep.ppm) > 20) {
					debugger;
				}
			}

			current_pep.PeptideID = pep.PeptideID;
			current_pep.FileID = pep.FileID;
			current_pep.Sequence = pep.Sequence;
			current_pep.ScoreValue = pep.ScoreValue;
			current_pep.mass = pep.Mass;
		}).then(function() {
			current_pep.modifications = [ [ 1, 1 , all_mods[current_pep.PeptideID] ] ];
			all_values[current_pep.FileID].push(derive_ppm_batch(current_pep));
			var foo = "";
			Object.keys(all_values).forEach(function(fileID) {
				foo = foo+all_values[fileID].map(function(pep) { return fileID+"\t"+pep.ppm+"\t"+pep.score; }).join("\n");
				foo = foo+"\n";
			});
			require('fs').writeFileSync('ppms.txt',foo);

			debugger;
		});
	});

};

exports.get_ppm = get_ppm;
exports.annotate_ppm = annotate_ppm;
exports.retrieve_all_ppms = retrieve_all_ppms;