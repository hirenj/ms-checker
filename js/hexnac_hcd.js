var util = require('util'),
    sqlite3 = require('sqlite3');
var spectra = require('./spectrum');
var peptide = require('./peptide');
var promisify_sqlite = require('./sqlite-promise');

var HexNAcHCD = function HexNAcHCD() {
};

util.inherits(HexNAcHCD,require('./processing-step.js'));

var open_db = function(filename) {
    return new Promise(function(resolve,reject) {
        var db = new sqlite3.Database(filename,sqlite3.OPEN_READONLY,function(err) {
            if (err) {
                reject(err);
                return;
            }
            resolve(db);
        });
    });
};

var get_db = function(filename) {
    return open_db(filename).then(function(db) {
        promisify_sqlite(db);
        db.partner = true;
        return db;
    });
};

module.exports = exports = new HexNAcHCD();

const default_error = {ppm: 15};

exports.error = default_error;

var check_mass = function(mass,mz,error) {
    if (! error ) {
        return mz == mass;
    }
    if (error.ppm) {
        return mz <= mass*(1 + error.ppm/1000000) && mz >= mass*(1 - error.ppm/1000000);
    }
    if (error.da) {
        return mz <= (mass+error.da) && mz >= (mass-error.da);
    }

    return false;
};

var is_galnac_mass = function(mz,error) {
    var galnac_138 = 138.055;
    var galnac_168 = 168.066;
    return (check_mass(galnac_138,mz,error) || check_mass(galnac_168,mz,error) );
};

var is_glcnac_mass = function(mz,error) {
    var glcnac_126 = 126.055;
    var glcnac_144 = 144.065;
    return (check_mass(glcnac_126,mz,error) || check_mass(glcnac_144,mz,error) );
};

var accept_mass = function(masses,intensities,mass,intensity) {
    var i = 0;
    var inserted = false;
    for (i = 0; ! inserted && i < masses.length; i++) {
        if (Math.abs(masses[i] - mass) <= 1) {
            if (intensity > intensities[i]) {
                intensities[i] = intensity;
                masses[i] = mass;
            }
            inserted = true;
        }
    }
    if ( ! inserted ) {
        masses.push(mass);
        intensities.push(intensity);
    }
};

var find_partner_hcd = function(spectrum,glyco_mass,db) {
    var self = this;

    if (typeof db === 'string') {
        return get_db(db).then(find_partner_hcd.bind(self,spectrum,glyco_mass));
    }

    if ( ! db.partner ) {
        return spectra.get_related_spectra(db,spectrum);
    }

    return spectra.match_spectrum_data(db, spectrum.scan, spectrum.rt, spectrum.charge, spectrum.mass - glyco_mass).then(function(spectra) {
        return spectra;
    });
};

var search_partner_hcd_in_msfs = function(spectrum,db,glyco_mass) {
    var self = this;
    // Search in db for another spectrum that matches this

    return Promise.all(
        [db].concat(self.sibling_msfs)
        .map(find_partner_hcd.bind(self,spectrum,glyco_mass))
    ).then(function(partners) {
        partners = Array.prototype.concat.apply([], partners).filter(function(spec) { return spec; });
        return partners;
    }).catch(function(err) {
        console.error(err);
    });
};

const mod_masses = {
    'HexNAc' : 203.079373,
    'HexHexNAc' : 365.132196,
    'Hex1HexNAc1' : 365.132196
};

var calculate_glyco_composition = function(composition,mods) {
    return composition.map(function(comp) {
        var bits = comp.split('x');
        var count = parseInt(bits.shift());
        var mass = mod_masses[bits.join('x')];
        if ( ! mass ) {
            console.log("Missing mass for ",bits.join('x'));
        }
        return count * mass;
    }).reduce(function(curr,next) { return curr+next; },0);
}

var decide_spectrum = function(db,pep,spectrum) {
    var self = this;
    if (! spectrum ) {
        return;
    }
    if (spectrum.activation !== 'HCD' && spectrum.activation !== 'CID' ) {
        var glyco_mass = calculate_glyco_composition(pep.Composition || peptide.composition(pep.modifications));
        return search_partner_hcd_in_msfs.bind(self)(spectrum,db,glyco_mass)
               .then(check_galnac_glcnac_ratio.bind(null,pep));
    }
    return check_galnac_glcnac_ratio(pep,spectrum);
};

var check_galnac_glcnac_ratio = function(pep,spectrum,partner) {
    if (! spectrum ) {
        return;
    }
    if (Array.isArray(spectrum)) {
        spectrum.forEach(function(spec) {
            check_galnac_glcnac_ratio(pep,spec,true);
        });
        return;
    }
    if (spectrum.activation !== 'HCD' && spectrum.activation !== 'CID' ) {
        return;
    }
    var galnac_intensities = [];
    var glcnac_intensities = [];
    var galnac_masses = [];
    var glcnac_masses = [];

    var error = this.error || {'ppm' : 15};
    if ( ! partner ) {
        pep.activation = spectrum.activation;
    }

    // Check to see what the ratio (mz-138 + mz-168) / (mz-126 + mz-144) is
    // if it is within 0.4 - 0.6, it is a GalNAc, 2.0 or greater, GlcNAc

    spectrum.peaks.forEach(function(peak) {
        var mass = peak.mass;
        var intensity = peak.intensity;
        if ( is_glcnac_mass(mass,error) ) {
            accept_mass( glcnac_masses, glcnac_intensities , mass, intensity);
        } else if ( is_galnac_mass(mass,error) ) {
            accept_mass( galnac_masses, galnac_intensities , mass, intensity);
        }
    });

    var galnac_count = galnac_masses.length;
    var glcnac_count = glcnac_masses.length;

    var galnac_intensity = galnac_intensities.reduce(function(curr,next) { return curr+next; },0);
    var glcnac_intensity = glcnac_intensities.reduce(function(curr,next) { return curr+next; },0);

    if (galnac_count < 2 || glcnac_count < 2 || galnac_count > 2 || glcnac_count > 2) {
        // console.log("Do not have 4 peaks to get the intensity for ",pep.SpectrumID,galnac_masses,glcnac_masses);
    }
    if (galnac_count == 2 && glcnac_count == 2) {
        var ratio = galnac_intensity / glcnac_intensity;
        var hexnac_type = (ratio <= 0.95) ? 'GalNAc' : ( (ratio >= 1.95) ? 'GlcNAc' : 'Unknown' );
        if (! pep.hexnac_type || pep.hexnac_type == 'Unknown') {
            pep.galnac_intensity = galnac_intensity;
            pep.glcnac_intensity = glcnac_intensity;
            pep.hexnac_type = hexnac_type;
        } else if (pep.hexnac_type && pep.hexnac_type != hexnac_type ) {
            console.log("Conflicting HexNAc types ", pep.hexnac_type,hexnac_type);
        }
    }
    return;
};

var guess_hexnac = function(db,sibling_msfs,peps) {
    var self = this;
    var to_cut = [].concat(peps);
    var split_length = 50;
    self.notify_task('Processing HexNAc HCD spectra');
    self.notify_progress(0,1);

    // We should only be procesing one set of peptides
    // at a time (i.e. one MSF file)
    self.sibling_msfs = sibling_msfs;

    var total = peps.length;
    result = Promise.resolve().then(function() {

        self.notify_progress(total-to_cut.length,total);

        if (to_cut.length < 1) {
            return peps;
        }
        return Promise.all( 
            to_cut.splice(0,split_length).map(function(pep) {
                return spectra.get_spectrum(db,pep).
                       then( decide_spectrum.bind(self,db,pep) );
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

var test_spectra = function(db,max_spectrum_id,error) {
    var self = this;
    self.error = error;
    var peps = [];
    while (max_spectrum_id > 0) {
        peps.push({'SpectrumID' : max_spectrum_id });
        max_spectrum_id -= 1;
    }
    return self.guess_hexnac(db,peps).then(function() {
        return peps.filter(function(pep) { return pep.activation; });
    }).then(function(filtered_peps) {
        console.log(filtered_peps);
        return filtered_peps;
    });
}

exports.guess_hexnac = guess_hexnac;
exports.test_spectra = test_spectra;