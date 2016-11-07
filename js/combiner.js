"use strict";

var peptide = require('./peptide');
var uniprot_meta = require('./uniprot').uniprot;
var nconf = require('nconf');
var quantitative = require('./quantitative');

var not_hcd_filter = function(pep) {
    return (pep.activation !== 'HCD' && pep.activation !== 'CID');
};

var onlyUnique = function(value, index, self) {
    return self.indexOf(value) === index;
};


var clone = function(objectToBeCloned) {
  // Basis.
  if (!(objectToBeCloned instanceof Object)) {
    return objectToBeCloned;
  }

  var objectClone;

  // Filter out special objects.
  var Constructor = objectToBeCloned.constructor;
  switch (Constructor) {
    // Implement other special objects here.
    case RegExp:
      objectClone = new Constructor(objectToBeCloned);
      break;
    case Date:
      objectClone = new Constructor(objectToBeCloned.getTime());
      break;
    default:
      objectClone = new Constructor();
  }

  // Clone each property.
  for (var prop in objectToBeCloned) {
    objectClone[prop] = clone(objectToBeCloned[prop]);
  }

  return objectClone;
};

var check_single_site = function(pep) {
    var indices = [];
    var str = pep.Sequence;
    for(var i=0; i<str.length;i++) {
        if (str[i].match(peptide.site_re)) indices.push(i+1);
    }
    if (pep.modifications || ! pep.Composition) {
        if ( ! pep.Composition && ! pep.modifications) {
            debugger;
        }
        return;
    }
    if (pep.Composition.length > 1) {
        return;
    }
    var total_sites = pep.Composition.map(function(comp) {  return parseInt(comp.split('x')[0]); }).reduce(function(old,n) { return old + n; },0);
    if (indices.length == total_sites && ! pep.modifications) {
        pep.modifications = indices.map(function(idx) {  return [ idx, pep.Composition[0].slice(2) ] ; });
        pep.made_ambiguous = 'inferred';
    }
};

var rewrite_deamidation = function(peps) {
    peps.forEach(function(pep) {
        (pep.modifications || []).forEach(function(mod) {
            if (mod[1].indexOf('Deamidated') >= 0) {
                mod[1] = 'GlcNAc(b1-4)GlcNAc';
            }
        })
    });
};

var combine_all_peptides = function(peps) {
    var data = {};
    var peptides_by_spectrum = {};
    peps.forEach(function(pep) {
        if ( ! peptides_by_spectrum[pep.SpectrumID]) {
            peptides_by_spectrum[pep.SpectrumID] = [];
        }
        peptides_by_spectrum[pep.SpectrumID].push(pep);
    });

    // Search by Spectrum ID:
    // Mark ambiguous peptides - if we have different sets of modifications for each spectrum, we can count that as ambiguous. Maybe mark out the possible amino acids? + total number of sites
    Object.keys(peptides_by_spectrum).forEach(function(spectrumID) {

        // Don't accept any peptides that has been identified using HCD
        // we don't trust the assignment, so we can skip them

        var peps = peptides_by_spectrum[spectrumID];

        // We want to count the number of modification keys for Non-HCD spectra
        // that identify sites. If there are multiple modification keys
        // it's a non-unique identification, so we should just use the
        // composition

        var mods = peps.filter(not_hcd_filter).filter( (pep) => ! pep.made_ambiguous ).map(function(pep) { return peptide.modification_key(pep); }).filter(onlyUnique).filter(function(key) { return key.indexOf('-') >= 0; });

        peps.forEach(function(pep) {
            if (mods.length > 1 || pep.activation == 'HCD' || pep.activation == 'CID') {
                if ( ! pep.Composition && pep.modifications && pep.modifications.length > 0 ) {
                    // We should be a bit smarter about this, grouping compositions appropriately
                    pep.Composition =  peptide.composition(pep.modifications);
                }
                if ( ! pep.modifications ) {
                    return;
                }
                if (mods.length > 1) {
                    console.log("We have multiple site assignment possibilities for a spectrum, and they haven't been removed by the fragmentation checker ",mods,peps,pep.modification_peptides);
                }
                pep.modifications.forEach(function(mod) {
                    mod[3] = 1;
                    mod[4] = pep.Sequence.length;
                });
                pep.possible_mods = pep.modifications;
                delete pep.modifications;
            }
        });

    });
    peptides_by_spectrum = null;

    var grouped_peptides = {};
    peps.forEach(function(pep) {

        check_single_site(pep);

        var pep_key = pep.Sequence+"-"+peptide.modification_key(pep);
        var current_uniprots = [];
        var current_starts = [];

        if ( ! grouped_peptides[pep_key] ) {
            grouped_peptides[pep_key] = [];
        } else {
            current_uniprots = grouped_peptides[pep_key][0].uniprot;
            current_starts = grouped_peptides[pep_key][0].starts || [grouped_peptides[pep_key][0].pep_start];
        }
        if (typeof pep.starts == 'undefined') {
            if (current_uniprots.indexOf(pep.uniprot[0]) < 0) {
                current_uniprots.push(pep.uniprot[0]);
                current_starts.push(pep.pep_start);
            }
        }
        for (var i = 0; i < (pep.starts || []).length; i++) {
            if (current_uniprots.indexOf(pep.uniprot[i]) < 0) {
                current_uniprots.push(pep.uniprot[i]);
                current_starts.push(pep.starts[i]);
            }
        }
        grouped_peptides[pep_key].push(pep);
        if (current_uniprots.length > 1) {
            grouped_peptides[pep_key][0].starts = current_starts;
        } else {
            grouped_peptides[pep_key][0].pep_start = current_starts[0];
        }
    });

    // Search by Peptide / (modifications + ambig key)
    // Combine quantitated peptides - seperating out into ambiguous / unambiguous.
    Object.keys(grouped_peptides).forEach(function(pep_key) {

        // HCD identification is not used for quantification

        var peps = grouped_peptides[pep_key].filter(not_hcd_filter);
        var quan_peps = peps.filter(function(pep) { return pep.acceptedquant && "QuanResultID" in pep;  });
        var singlet_peps = quan_peps.filter(function(pep) { return pep.QuanChannelID.length < 2 && pep.has_pair == false; });
        var paired_singlets = quan_peps.filter(function(pep) { return pep.QuanChannelID.length < 2 && pep.has_pair == true; });
        var ratioed_peps = quan_peps.filter(function(pep) { return pep.QuanChannelID.length == 2; });
        var target_ratio = null;
        var target_ratio_mad = null;
        if (ratioed_peps.length > 0) {
            var seen_quan_result_ids = {};

            // We wish to consider each ratio once for each
            // quantitation result. If there are multiple modifications
            // or some other thing that boosts the number of peptides
            // this will ignore them.
            //
            // There are cases where there is a single QuanResultID,
            // but there are multiple spectra (and identified peptides) associated
            // with this one QuanResultID

            var all_ratios = ratioed_peps.filter(function(pep) {
                var seen = pep.QuanResultID in seen_quan_result_ids;
                seen_quan_result_ids[pep.QuanResultID] = 1;
                return ! seen;
            }).map(function(pep) { return pep.areas[quantitative.ko_channel] / pep.areas[quantitative.wt_channel]; });

            target_ratio = Math.median( all_ratios );

            // We wish to record the mad to see if we are picking up
            // different populations of sites. Higher mad means that
            // there that there can be another population of values
            // here. This is especially useful for the ambiguous
            // modification sites, where we might see various populations.

            // We can always go back to them to examine them later.

            target_ratio_mad = Math.median(all_ratios.map(function(ratio) { return Math.abs(ratio - target_ratio); }));


            singlet_peps.forEach(function(pep) { pep.used = false; });
            paired_singlets.forEach(function(pep) { pep.used = false; });
            ratioed_peps.forEach(function(pep) { pep.used = true; });

        } else if (singlet_peps.length > 0)  {

            target_ratio = singlet_peps[0].QuanChannelID[0] == quantitative.wt_channel ? 1/100000 : 100000;

            var channel_ids = singlet_peps.map(function(pep) { return pep.QuanChannelID[0]; }).filter(onlyUnique);
            singlet_peps.forEach(function(pep) {
                pep.used = true;
            });
            if (channel_ids.length > 1) {
                target_ratio = 'conflicting_singlets';
            }
        } else {
            paired_singlets.forEach(function(pep) {
                pep.used = true;
            });
        }
        if (target_ratio) {
            quan_peps.forEach(function(pep) {
                if (target_ratio_mad !== null) {
                    pep.CalculatedRatio_mad = target_ratio_mad;
                }
                pep.CalculatedRatio = target_ratio;
            });
        }

        peps.filter(function(pep) { return (! pep.acceptedquant) && "QuanResultID" in pep;  }).forEach(function(pep) {
            // For peptides that aren't accepted as being quantitated by Proteome Discoverer
            // we should only report them if they don't exist in a group along side
            // any other peptides that have a ratio.
            if (! target_ratio) {
                if (pep.QuanChannelID.length > 1) {
                    pep.unaccepted_quant = "ratio";//pep.areas[quantitative.ko_channel] / pep.areas[quantitative.wt_channel];
                } else if (pep.QuanChannelID.length == 1) {
                    pep.unaccepted_quant = "singlet_"+pep.QuanChannelID[0];
                }
            }
        });


    });

    peps.forEach(function(pep) {
        pep.uniprot = uniprot_meta.updateIds(pep.uniprot);
        peptide.fix_site_numbers(pep);
    });


    Object.keys(grouped_peptides).forEach(function(pep_key) {
        var peps = grouped_peptides[pep_key];
        if (peps.length < 1) {
            return;
        }
        var first_pep = peps[0];
        var quant = null;
        var high_sn = false;
        var confident_singlet = false;
        var has_possible_mods = false;
        var is_singlet = false;
        var hexnac_type = {};
        var hexnac_ratios = [];
        var max_score = null;
        var unaccepted_quant = null;


        var spectra = [];
        var activations = [];
        var areas = {'medium': [], 'light' : []};
        var quant_intensity = {'min' : null, 'max' : null};
        peps.forEach(function(pep) {
            if (pep.acceptedquant) {
                if ("CalculatedRatio" in pep) {
                    quant = pep.CalculatedRatio;
                    high_sn = true;
                    if (pep.CalculatedRatio == 1/100000 || pep.CalculatedRatio == 100000) {
                        is_singlet = true;
                    }
                } else {
                    is_singlet = true;
                }
            }
            if ("has_low_sn" in pep && ! pep.has_low_sn) {
                high_sn = true;
            }
            if ("has_high_sn" in pep && pep.has_high_sn) {
                confident_singlet = true;
            }
            if (pep.acceptedquant && "has_pair" in pep && pep.has_pair === true && pep.used && ( pep.activation !== 'HCD' && pep.activation !== 'CID' )) {
                // Potential ratio 1/100000 in the potential_light
                // Potential ratio 100000 in the potential_medium 
                quant = pep.QuanChannelID[0] == quantitative.wt_channel ? 'potential_wt' : 'potential_ko';
            }

            if ( ! pep.acceptedquant ) {
                unaccepted_quant = pep.unaccepted_quant;
            }

            if ("hexnac_type" in pep) {
                hexnac_type[pep.hexnac_type] = true;
                if (("GlcNAc" in hexnac_type || "GalNAc" in hexnac_type)) {
                    delete hexnac_type['unclear'];
                    delete hexnac_type['ND'];
                }
                if (pep.hexnac_type == 'ND') {
                    hexnac_ratios.push('ND');
                } else {
                    hexnac_ratios.push((pep.galnac_intensity/pep.glcnac_intensity).toFixed(2));
                }
            }
            if ("possible_mods" in pep && pep.activation !== 'HCD') {
                has_possible_mods = true;
            }

            // The areas arrays should only contain the areas for the
            // quantitations that have been accepted.

            if (pep.areas && pep.acceptedquant && pep.areas['medium']) {
                areas['medium'].push(pep.areas['medium']);
                if ( ! pep.areas['light']) {
                    areas['light'].push(0);
                }
            }
            if (pep.areas && pep.acceptedquant && pep.areas['light']) {
                areas['light'].push(pep.areas['light']);
                if ( ! pep.areas['medium']) {
                    areas['medium'].push(0);
                }
            }
            if ( ! pep.areas || ! (pep.areas['light'] || pep.areas['medium']) || ! pep.acceptedquant ) {
                areas['light'].push(0);
                areas['medium'].push(0);
            }
            if (pep.quant_intensity_range) {
                if (quant_intensity.min === null || pep.quant_intensity_range.min < quant_intensity[0]) {
                    quant_intensity.min = pep.quant_intensity_range[0];
                }
                if (quant_intensity.max === null || pep.quant_intensity_range.max > quant_intensity[1]) {
                    quant_intensity.max = pep.quant_intensity_range[1];
                }
            }
            spectra.push( {'score' : parseFloat(pep.score.toFixed(2)), 'rt' : parseFloat(pep.retentionTime.toFixed(2)), 'scan' : pep.scan, 'ppm' : parseFloat(pep.ppm.toFixed(2)), 'mass' : parseFloat(pep.mass.toFixed(2)) , 'charge' : pep.charge } );
            activations.push( pep.activation );
        });

        // If we have a ratio / singlet and an unaccepted quant
        // we should just accept the quant, and ignore the unaccepted quant info
        // There is the possibility of checking to make sure that the rejected
        // singlet goes in the same direction as the other quants, but
        // I think the logic of that is a bit iffy.

        if (quant && unaccepted_quant) {
            unaccepted_quant = null;
        }

        var block = {
            'sequence' : first_pep.Sequence
        };

        if (first_pep.uniprot.length > 1) {
            block.multi_protein = "protein";
        }
        if (first_pep.multi_peptide) {
            block.multi_protein = block.multi_protein ? block.multi_protein+',peptide' : 'peptide';
        }

        if (quant !== null) {
            block.quant = isNaN(parseInt(quant)) ? {'quant' : quant } : {'quant' : quant, 'mad' : first_pep.CalculatedRatio_mad };
            block.quant.channels = quantitative.ko_channel[0].toUpperCase()+':'+quantitative.wt_channel[0].toUpperCase();
            if (is_singlet) {
                block.quant.singlet_confidence = confident_singlet ? 'high' : 'low';
                if ( ! high_sn ) {
                    block.quant.singlet_confidence = 'low_sn';
                }
                block.quant.intensities = quant_intensity;
            }
        }
        if (unaccepted_quant) {
            block.quant = {'quant' : "rejected_"+unaccepted_quant };
        }
        if (Object.keys(hexnac_type).length > 0) {
            block.hexnac_type = Object.keys(hexnac_type).sort();
            block.hexnac_ratio = hexnac_ratios.join(',');
        }

        if (nconf.get('feature_enable_nglyco')) {
            rewrite_deamidation([first_pep]);
        }

        if (first_pep.modifications) {
            block.sites = first_pep.modifications.map(function(mod) { return [ mod[0], mod[1] ]; } );
        }

        if (first_pep.Composition) {
            block.composition = first_pep.Composition;
        } else if (first_pep.modifications) {
            var count = 0;
            block.composition = peptide.composition(first_pep.modifications);
        }

        block.peptide_start = first_pep.pep_start + 1;

        block.spectra = spectra;
        block.activation = activations.filter(onlyUnique);
        if (block.quant) {
            block.quant_areas = areas;
        }

        if (has_possible_mods) {
            block.ambiguous_mods = peps.filter(function(pep) { return pep.possible_mods; }).map(function(pep) {
                var written_mods = write_possible_mods(pep.possible_mods);
                var temp_block = {'ambiguous_mods' : [written_mods] };
                rescale_ids(pep.pep_start, first_pep.pep_start, temp_block);
                return temp_block.ambiguous_mods[0];
            }).filter(onlyUnique);
            block.sites = block.ambiguous_mods[0].split(',').filter((mod) => mod.indexOf('-') < 0).map((mod) => { let mod_bits = mod.split(/[\(\)]/).slice(0,2); mod_bits[0] = +mod_bits[0]; return mod_bits; });
            block.ambiguous_sites = consolidate_ambiguous_compositions(block.ambiguous_mods[0].split(',').filter((mod) => mod.indexOf('-') > 0));

            delete block.ambiguous_mods;


            block.made_ambiguous = '';

            var ambig_reasons = peps.filter(function(pep) { return "made_ambiguous" in pep; }).map(function(pep) { return pep.made_ambiguous == true ? 'missing_site_coverage' : pep.made_ambiguous; });
            if (ambig_reasons.length > 0) {
                block.made_ambiguous = ambig_reasons.sort().filter(onlyUnique).join(',');
            }
        } else {
            block.made_ambiguous = peps.filter(function(pep) { return pep.made_ambiguous; }).map(function(pep) { return pep.made_ambiguous; }).filter(onlyUnique)[0];
            if (! block.made_ambiguous) {
                delete block.made_ambiguous;
            }
        }
        var genes = first_pep.gene || [];

        first_pep.uniprot.forEach(function(uniprot,idx) {
            if ( ! data[uniprot] ) {
                data[uniprot] = [];
            }
            var block_copy = clone(block);
            if (first_pep.starts && first_pep.pep_start != first_pep.starts[idx]) {
                rescale_ids(first_pep.pep_start,first_pep.starts[idx],block_copy);
            }
            if (genes.length) {
                block_copy.gene = genes.shift();
            }
            data[uniprot].push(block_copy);
        });
    });
    return data;
};

var rescale_ids = function(old,new_pos,peptide) {
    if (old == new_pos) {
        return;
    }
    var ambiguous = peptide.ambiguous_mods;
    if (ambiguous) {
        ambiguous = ambiguous.map(function(ambig_pep) {
            return ambig_pep.split(',').map(function(ambig) {
                var bits = ambig.split('(');
                var positions = bits[0].split('-');
                var mapped = positions.map(function(pos) {
                    return parseInt(pos) - old + new_pos;
                });
                return mapped.join('-')+"("+bits[1];
            }).join(',');
        });
        peptide.ambiguous_mods = ambiguous;
    }
    if (peptide.sites) {
        peptide.sites.forEach(function(comp) {
            comp[0] = parseInt(comp[0]) - old + new_pos;
        });
    }
    peptide.peptide_start = new_pos + 1;
}

var consolidate_ambiguous_compositions = function(sites) {
    let site_counts = {};
    let all_sites = [];
    sites.forEach( site => {
        site_counts[site] = (site_counts[site] || 0) + 1;
        all_sites.push(site);
    });
    return all_sites.filter(onlyUnique).map( mod => {
        let mod_bits = mod.split(/[\(\)]/).slice(0,2);
        let range_ends = mod_bits[0].split('-').map( end => +end );
        if (site_counts[mod] > 1) {
            mod_bits[1] = site_counts[mod] + 'x' + mod_bits[1];
        }
        return [ range_ends, mod_bits[1] ];
    });
};

var write_possible_mods = function(mods) {
    var mod_compositions = {};
    return mods.sort(function(a,b) {
        var sort_val_a = a[3] ? a[3] : a[0];
        var sort_val_b = b[3] ? b[3] : b[0];
        if (sort_val_a == sort_val_b) {
            sort_val_a = a[4] ? a[4] : a[0];
            sort_val_b = b[4] ? b[4] : b[0];
        }
        return sort_val_a - sort_val_b;
    }).map(function(mod) {
        if (mod[3]) {
            if (mod[3] == mod[4]) {
                return ""+mod[3]+"("+mod[1]+")";
            }
            return ""+mod[3]+"-"+mod[4]+"("+mod[1]+")";
        }
        return mod[0]+"("+mod[1]+")";
    }).join(',');
};


exports.combine = combine_all_peptides;