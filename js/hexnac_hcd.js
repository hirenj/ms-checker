var util = require('util');
var spectra = require('./spectrum');

var HexNAcHCD = function HexNAcHCD() {
};

util.inherits(HexNAcHCD,require('./processing-step.js'));

module.exports = exports = new HexNAcHCD();


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
    if (spectrum.activation !== 'HCD' && spectrum.activation !== 'CID' ) {
        return;
    }
    pep.activation = spectrum.activation;

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

var guess_hexnac = function(db,peps) {
    var self = this;
    var to_cut = [].concat(peps);
    var processing_node = null;
    var result = spectra.init_spectrum_processing_num(db).then(function(node) {
        processing_node = node;
        if (typeof self.conf.get('hcd-processing-node') !== 'undefined') {
            processing_node = parseInt(self.conf.get('hcd-processing-node'));
        }
        if ( processing_node !== 0 && ! processing_node ) {
            throw new Error("No Processing node for HCD");
        }
        console.log("We only want spectra from Processing Node ",processing_node ? processing_node : 'Any processing node');
        return true;
    });
    var split_length = 50;

    self.notify_task('Processing HexNAc HCD spectra');
    self.notify_progress(0,1);

    var total = peps.length;
    result = result.then(function() {

        self.notify_progress(total-to_cut.length,total);

        if (to_cut.length < 1) {
            return peps;
        }
        return Promise.all(  to_cut.splice(0,split_length).map(function(pep) {
                                return spectra.get_spectrum(db,pep,processing_node).then( check_galnac_glcnac_ratio.bind(null,pep) );
                            }) ).then(arguments.callee);
    }).catch(function(err) {
        if (err.message === "No Processing node for HCD") {
            self.notify_progress(1,1);
            return Promise.resolve(peps);
        }
        throw err;
    });
    return result;
};

exports.guess_hexnac = guess_hexnac;