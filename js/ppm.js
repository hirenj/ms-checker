var masses = require('./fragmentions').calculate_peptide_mass;
var quantitative = require('./quantitative');

// We love Protein Prospector when trying to
// figure out what the right masses should be
// http://prospector.ucsf.edu/prospector/cgi-bin/mssearch.cgi


var calculate_theoretical_mass = function(peptide) {
	var aa_masses = masses(peptide);
	// We want to take the Mass column from SpectrumHeaders (which is the MH+) in the spreadsheets
	peptide.mass_theoretical = (aa_masses.reduce(function(p,n) { return p+n; },0));
};

var get_ppm = function(peptide) {
	calculate_theoretical_mass(peptide);
	peptide.ppm  = 1e06*((peptide.mass / peptide.mass_theoretical) - 1);
}

var annotate_ppm = function(peps) {
	peps.forEach(get_ppm);
	return peps;
}

exports.get_ppm = get_ppm;
exports.annotate_ppm = annotate_ppm;