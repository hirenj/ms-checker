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
WHERE SpectrumHeaders.Charge BETWEEN 2 AND 3 \
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
		// console.log(3246115 - mods.PeptideID);
		if (mods.PeptideID !== last_peptide_id) {
			all_mods[last_peptide_id] = all_mods[last_peptide_id].reduce(function(p,n) { return p + n[2] },0);
			last_peptide_id = mods.PeptideID;
		}
	};

	var all_mods_promise = db.each(retrieve_mods_sql,[],mods_function);

	//ggplot(foo)+geom_boxplot(aes(x=cuts,y=V3),coef=3)
	//foo$cuts = cut(foo$V2, breaks = seq(-15, 15, by = 0.5))
	//blah = by(foo$V3,cut(foo$V2,breaks=seq(-15,15,by=0.5)),function(x) {length(x[x > 3.104]) })
	//blah = by(foo$V3,cut(foo$V2,breaks=seq(-15,15,by=0.5)),function(x) { browser(); length(x[x > quantile(x,0.99)]) })

	return all_mods_promise.then(function() {
		return db.each(retrieve_ppm_all_peptides_sql,[],function(err,pep) {
			if ( ! all_values[pep.FileID]) {
				all_values[pep.FileID] = [];
			}
 			// console.log(3246115 - pep.PeptideID);
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

			// var foo = "";
			// foo = foo+all_values[fileID].map(function(pep) { return fileID+"\t"+pep.ppm+"\t"+pep.score; }).join("\n");
			// foo = foo+"\n";
			// // Plot values
			// require('fs').writeFileSync('ppms.txt',foo);

			var total_values = [];
			Object.keys(all_values).forEach(function(fileID) {
				total_values = total_values.concat(all_values[fileID]);
			});

			// Bin the scores by 0.5 ppm windows

			var bin_info = bin_values(total_values);

			var bins = bin_info.bins;

			var first_bin = bin_info.first;

			// Extract the count of outliers in each bin in top 99.9%

			var quantiles = extract_quantiles(bins);

			var total_sum = quantiles.reduce(function(p,n) { return p+n; },0);

			// Look for window where the number of outliers is > 9/10*(total sum/100)

			var left_bin = first_bin;

			while (quantiles[0] !== null && quantiles.length > 2 && (quantiles[0] + quantiles[1] + quantiles[2])/3 < 0.9*(total_sum/100)) {
				left_bin += 0.5;
				quantiles.shift();
			}

			quantiles = quantiles.reverse();
			var right_bin = first_bin + bins.length*0.5;
			while (quantiles[0] !== null && quantiles.length > 2 && (quantiles[0] + quantiles[1] + quantiles[2])/3 < 0.9*(total_sum/100)) {
				right_bin -= 0.5;
				quantiles.shift();
			}

		});
	});

};

var bin_values = function(values) {
	var ppms = values.map(function(point) { return point.ppm; }).sort(function(a,b) { return a-b; });
	var min_ppm = ppms[0];
	min_ppm = Math.floor(2*min_ppm)/2;
	var bins = [];
	values.forEach(function(point) {
		var bin_idx = Math.floor(2*(point.ppm - min_ppm));
		if ( ! bins[bin_idx] ) {
			bins[bin_idx] = [];
		}
		bins[bin_idx].push(point.score);
	});
	return { bins: bins, first: min_ppm };
}

var quantile = function(array,p) {
	var sorted = array.slice().sort(function(a,b) { return a - b});
	var idx = sorted.length * p;
	return sorted[Math.ceil(idx) - 1];
}

var extract_quantiles = function(bins) {
	var max_value = 0.5*(quantile(bins[1],0.999)+quantile(bins[2],0.999));
	console.log(max_value);
	return bins.map(function(bin) { return bin.filter(function(score) { return score > max_value; }).length });
};

var cumulative_sum = function(bins) {
	var total = 0;
	return bins.map(function(count) { total = total + count; return total; });
};

var diff = function(array) {
	var result = [];
	for (var i = 0; i < (array.length - 1); i++) {
		result[i] = array[i+1] - array[i];
	}
	return result;
}


exports.get_ppm = get_ppm;
exports.annotate_ppm = annotate_ppm;
exports.retrieve_all_ppms = retrieve_all_ppms;