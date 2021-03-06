var util = require('util');
var spectra = require('./spectrum');
var quantitative = require('./quantitative');
var peptide = require('./peptide');

var site_re = peptide.site_re;

const MASS_H = 1.007825035;
const MASS_ELECTRON = 0.00054858026;
const MASS_H_MINUS_ELECTRON = 1.007276455;
const MASS_C = 12.000000;
const MASS_N = 14.003074;
const MASS_O = 15.99491463;

const MASS_AMINO_ACIDS = {
    A:71.03712,
    B:114.53494,
    C:103.00919,
    D:115.02695,
    E:129.0426,
    F:147.06842,
    G:57.02147,
    H:137.05891,
    I:113.08407,
    K:128.09497,
    L:113.08407,
    M:131.0405,
    N:114.04293,
    P:97.05277,
    Q:128.05858,
    R:156.10112,
    S:87.03203,
    T:101.04768,
    V:99.06842,
    W:186.07932,
    X:0.0,
    Y:163.06332,
    Z:128.55059
};

const retrieve_matching_node_config_sql = 'SELECT \
    ProcessingNodeNumber, ParameterName, ParameterValue \
FROM Peptides \
JOIN \
    ProcessingNodeParameters USING(ProcessingNodeNumber) \
WHERE \
    PeptideID = ? \
AND \
    Category = "4. Spectrum Matching" \
';

var FragmentIons = function FragmentIons() {
};

util.inherits(FragmentIons,require('./processing-step.js'));

module.exports = exports = new FragmentIons();

var onlyUnique = function(value, index, self) {
    return self.indexOf(value) === index;
};

var cached_node_filters = {};

var get_ion_matching_config = function(db,pep) {
    return db.all(retrieve_matching_node_config_sql,[pep.PeptideID]).then(function(configs) {
        var filters = [];
        if (configs.length < 1) {
            return function() { return false; };
        }
        if (cached_node_filters[configs[0].ProcessingNodeNumber]) {
            return cached_node_filters[configs[0].ProcessingNodeNumber];
        }
        configs.forEach(function(config) {
            if (config.ParameterValue == "1" || config.ParameterValue == "True" ) {
                return;
            }
            if (config.ParameterName == 'UseNeutralAIons') {
                filters.push(/a_[nh]/);
            }
            if (config.ParameterName == 'UseNeutralBIons') {
                filters.push(/b_[nh]/);
            }
            if (config.ParameterName == 'UseNeutralYIons') {
                filters.push(/y_[nh]/);
            }
            if (config.ParameterName == 'IonSerieA') {
                filters.push(/a/);
            }
            if (config.ParameterName == 'IonSerieB') {
                filters.push(/b/);
            }
            if (config.ParameterName == 'IonSerieC') {
                filters.push(/c/);
            }
            if (config.ParameterName == 'IonSerieX') {
                filters.push(/x/);
            }
            if (config.ParameterName == 'IonSerieY') {
                filters.push(/y/);
            }
            if (config.ParameterName == 'IonSerieZ') {
                filters.push(/z/);
            }
        });
        cached_node_filters[configs[0].ProcessingNodeNumber] = function(type) {
            return filters.reduce(function(val,regex) { return val ? val : regex.test(type); } , false );
        };
        return cached_node_filters[configs[0].ProcessingNodeNumber];
    });
};

var assign_peptide_ions = function(db,pep,debug) {
    return Promise.all( [ spectra.get_spectrum(db,pep), get_ion_matching_config(db,pep) ]).then(function(spec_data) {
        var spectrum = spec_data[0];
        var ion_filters = spec_data[1];
        if (debug) {
            debugger;
        }
        if ( ! spectrum ) {
            return [];
        }

        var theoretical_ions = calculate_fragment_ions(pep,spectrum.charge || 1).filter(function(ion) {
            return ! ion_filters(ion.type);
        });
        theoretical_ions.forEach(function(ion) {
            var min_mz = ion.mz - 0.01;
            var max_mz = ion.mz + 0.01;
            var matching_peaks = spectrum.peaks.filter(function(peak) {  return peak.mass <= max_mz && peak.mass >= min_mz;   });
            ion.peaks = matching_peaks;
        });
        return theoretical_ions.filter(function(ion) { return (ion.peaks || []).length > 0; });
    });
};

var filter_assigned_for_isotope_envelope = function(ions) {
    return ions.filter(function(ion) {
        var peaks = (ion.peaks || []).filter(function(peak) { return peak.charge !== 0 && peak.charge == ion.z; });
        if (peaks.length > 0) {
            return true;
        }
        return false;
    });
};

var get_b_ion_coverage = function(db,pep) {
    return assign_peptide_ions(db,pep).then(filter_assigned_for_isotope_envelope).then(function(ions) {

        var b_positions = ions.map( function(ion_data) {
            var ion = ion_data.type;
            var match;
            if (match = ion.match(/[xyz](?:_.+_)?(\d+)$/)) {
                return pep.Sequence.length - parseInt(match[1]);
            }
            if (match = ion.match(/[acb](?:_.+_)?(\d+)/)) {
                return parseInt(match[1]);
            }
        });

        return b_positions.filter(function(pos) { return pos && pos > 0; }).filter(onlyUnique).sort(function(a,b) { return a-b; });
    });
};

var theoretical_ions = function(db,pep,charge) {
    return get_ion_matching_config(db,pep).then(function(ion_filter) {
        return calculate_fragment_ions(pep,charge || 1).filter(function(ion) {
            return ! ion_filter(ion.type);
        });
    });
};

var matched_ions = function(db,pep,debug) {
    return assign_peptide_ions(db,pep,debug).then(filter_assigned_for_isotope_envelope);
};

var get_coverage_for_sites = function(pep,b_ions) {
    if (! pep.modifications || pep.modifications.length == 0 || b_ions.length == 0) {
        return [];
    }
    return pep.modifications.map(function(mod) {
       var pos = mod[0];
       var min = [0].concat(b_ions).filter(function(ion) { return ion < pos; }).reverse()[0];
       var max = b_ions.concat([pep.Sequence.length]).filter(function(ion) { return ion >= pos; })[0];
       return ""+min+"-"+pep.Sequence.substring(min,max);
    });
};

var resolve_glyco_modifications = function(pep) {
    return resolve_modifications(pep,site_re);
};

var resolve_modifications = function(pep,aa_re) {
    var mods = [].concat(pep.modifications);
    var groups = {};
    pep.modification_peptides.forEach(function(aas) {
        if ( ! groups[aas]) {
            groups[aas] = [];
        }
        groups[aas].push(mods.shift());
    });

    var is_ambiguous = false;
    Object.keys(groups).forEach(function(aagroup) {
        var seq = aagroup.split('-')[1];
        var start_pos = parseInt(aagroup.split('-')[0]);
        var start = null;
        var end;
        var match;
        if ( ! seq.match(aa_re) ) {
            made_ambiguous = true;
            pep.invalid_modification = true;
            groups[aagroup].forEach(function(mod) {
                mod[3] = 0;
                mod[4] = 0;
            });
            return pep;
        }
        if ( seq.match(aa_re).length > groups[aagroup].length ) {
            is_ambiguous = true;
            while ( match = aa_re.exec(seq) ) {
                if ( start === null ) {
                    start = start_pos + match.index;
                } else {
                    end = start_pos + match.index;
                }
            }
            groups[aagroup].forEach(function(mod) {
                mod[3] = start + 1;
                mod[4] = end + 1;
            });
        }
    });
    if (is_ambiguous) {
        if ( ! pep.Composition && is_ambiguous ) {
            pep.Composition =  peptide.composition(pep.modifications);
        }
        pep.possible_mods = pep.modifications;
        delete pep.modifications;
        pep.made_ambiguous = true;
    }
    return pep;
};

var batch_promise = function(list,split,mapper) {
    exports.notify_progress(0,1);
    var total = list.length;
    var to_cut = [].concat(list);
    var result = Promise.resolve(true);

    result = result.then(function() {
        exports.notify_progress( (total - to_cut.length), total );

        if (to_cut.length < 1) {
            return list;
        }
        return Promise.all( to_cut.splice(0,split).map(mapper) ).then(arguments.callee);
    });

    return result;
};

var validate_peptide_coverage = function(db,peptides) {
    var self = this;
    var total = null;
    console.time('Fragmentation');
    exports.notify_task('Validating fragmentation spectra');

    var wanted_peptides = peptides.filter(function(pep) { return ((pep.activation || "") !== "HCD") &&   ((pep.activation || "") !== "CID") && pep.modifications && pep.modifications.length > 0; });

    return batch_promise(wanted_peptides,50,function(pep,idx) {
        return get_b_ion_coverage(db,pep).then(function(ions) {
            pep.modification_peptides = get_coverage_for_sites(pep,ions);
            return pep;
        });
    }).then(function(modified_peps) {
        console.timeEnd('Fragmentation');
        modified_peps.forEach(resolve_glyco_modifications);
        exports.notify_progress('progress',1,1);
    }).then(function() {
        return peptides.filter(function(pep) { return ! pep.invalid_modification; });
    });
};

var sum_mods = function(array,idx) {
    if (array.length == 0) {
        return 0;
    }
    // Should we have unique masses here?
    return array.filter(function(mod) { return mod[0] === idx; }).map(function(mod) { return mod[2]; }).reduce(function(curr,next) { return curr + next; },0);
};

var calculate_peptide_masses = function(pep,modifications) {
    var mods = modifications || quantitative.modifications_cache[pep.PeptideID] || [];
    var masses = pep.Sequence.split('').map(function(aa,idx) { return MASS_AMINO_ACIDS[aa] + sum_mods(mods,idx+1);  });
    return masses.concat(MASS_H_MINUS_ELECTRON + 2*MASS_H + MASS_O);
};

// We love Protein Prospector when trying to
// figure out what the right masses should be
// http://prospector.ucsf.edu/prospector/cgi-bin/mssearch.cgi

var calculate_fragment_ions = function(pep,spectrum_charge) {
    var mods = quantitative.modifications_cache[pep.PeptideID] || [];
    var masses = pep.Sequence.split('').map(function(aa,idx) { return MASS_AMINO_ACIDS[aa] + sum_mods(mods,idx+1);  });
    var ions = [];
    var charge = spectrum_charge;
    var b_ion_base, y_ion_base;
    while(charge > 0) {
        b_ion_base = 0;

        // Should be minus electron here?
        y_ion_base = 2*MASS_H + MASS_O;

        for (var i = 0 ; i < pep.Sequence.length; i++) {
            b_ion_base += masses[i];
            y_ion_base += masses[(masses.length - 1) - i];


            ions = ions.concat( [
            { 'type' : 'b'+(i+1), 'mz' : (b_ion_base + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'b_h2o_'+(i+1), 'mz' : (b_ion_base - MASS_O - 2 * MASS_H + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'b_nh3_'+(i+1), 'mz' : (b_ion_base - MASS_O - MASS_N - 3* MASS_H + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'a'+(i+1), 'mz' : (b_ion_base - MASS_O - MASS_C + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'a_nh3_'+(i+1), 'mz' : (b_ion_base - MASS_O - MASS_C - MASS_N - 3 * MASS_H + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'a_h20_'+(i+1), 'mz' : (b_ion_base - 2*MASS_O - MASS_C - 2 * MASS_H + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'c'+(i+1), 'mz' : (b_ion_base + MASS_N + 3 * MASS_H + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'c_+1_'+(i+1), 'mz' : (b_ion_base + MASS_N + 3 * MASS_H + MASS_H + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'c_-1_'+(i+1), 'mz' : (b_ion_base + MASS_N + 3 * MASS_H - MASS_H + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'c_+2_'+(i+1), 'mz' : (b_ion_base + MASS_N + 3 * MASS_H + 2*MASS_H + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'c_-2_'+(i+1), 'mz' : (b_ion_base + MASS_N + 3 * MASS_H - 2*MASS_H + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'c_+3_'+(i+1), 'mz' : (b_ion_base + MASS_N + 3 * MASS_H + 3*MASS_H + charge * MASS_H) / charge, 'z' : charge },
            { 'type' : 'c_-3_'+(i+1), 'mz' : (b_ion_base + MASS_N + 3 * MASS_H - 3*MASS_H + charge * MASS_H) / charge, 'z' : charge },


            { 'type' : 'y'+(i+1), 'mz' : ((y_ion_base + charge * MASS_H) / charge) - MASS_ELECTRON, 'z' : charge },
            { 'type' : 'y_nh3_'+(i+1), 'mz' : ((y_ion_base - MASS_N - 3 * MASS_H + charge * MASS_H) / charge) - MASS_ELECTRON, 'z' : charge },
            { 'type' : 'y_h2o_'+(i+1), 'mz' : ((y_ion_base - 2* MASS_H - MASS_O + charge * MASS_H) / charge) - MASS_ELECTRON, 'z' : charge },
            { 'type' : 'x'+(i+1), 'mz' : ((y_ion_base + MASS_C + MASS_O - 2 * MASS_H + charge * MASS_H) / charge) - MASS_ELECTRON, 'z' : charge },
            { 'type' : 'z'+(i+1), 'mz' : ((y_ion_base - MASS_N - 2 * MASS_H + charge * MASS_H) / charge) - MASS_ELECTRON, 'z' : charge },
            { 'type' : 'z_+1_'+(i+1), 'mz' : ((y_ion_base - MASS_N - 2 * MASS_H + MASS_H + charge * MASS_H) / charge) - MASS_ELECTRON, 'z' : charge },
            { 'type' : 'z_-1_'+(i+1), 'mz' : ((y_ion_base - MASS_N - 2 * MASS_H - MASS_H + charge * MASS_H) / charge) - MASS_ELECTRON, 'z' : charge },
            { 'type' : 'z_+2_'+(i+1), 'mz' : ((y_ion_base - MASS_N - 2 * MASS_H + 2*MASS_H + charge * MASS_H) / charge) - MASS_ELECTRON, 'z' : charge },
            { 'type' : 'z_-2_'+(i+1), 'mz' : ((y_ion_base - MASS_N - 2 * MASS_H - 2*MASS_H + charge * MASS_H) / charge) - MASS_ELECTRON, 'z' : charge },
            { 'type' : 'z_+3_'+(i+1), 'mz' : ((y_ion_base - MASS_N - 2 * MASS_H + 3*MASS_H + charge * MASS_H) / charge) - MASS_ELECTRON, 'z' : charge },
            { 'type' : 'z_-3_'+(i+1), 'mz' : ((y_ion_base - MASS_N - 2 * MASS_H - 3*MASS_H + charge * MASS_H) / charge) - MASS_ELECTRON, 'z' : charge }

            ]);
        }
        charge -= 1;
    }
    return ions;
};
exports.clear_caches = function() { cached_node_filters = {}; };
exports.validate_peptide_coverage = validate_peptide_coverage;
exports.matched_ions = matched_ions;
exports.theoretical_ions = theoretical_ions;
exports.all_theoretical_ions = calculate_fragment_ions;
exports.calculate_peptide_mass = calculate_peptide_masses;


