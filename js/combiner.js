var peptide = require('./peptide');
var uniprot_meta = require('./uniprot');

var not_hcd_filter = function(pep) {
    return (pep.activation !== 'HCD' && pep.activation !== 'CID');
};

var onlyUnique = function(value, index, self) {
    return self.indexOf(value) === index;
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

        var mods = peps.filter(not_hcd_filter).map(function(pep) { return peptide.modification_key(pep); }).filter(onlyUnique).filter(function(key) { return key.indexOf('-') >= 0; });

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
        var pep_key = pep.Sequence+"-"+peptide.modification_key(pep);
        if ( ! grouped_peptides[pep_key] ) {
            grouped_peptides[pep_key] = [];
        }
        grouped_peptides[pep_key].push(pep);
    });

    // Search by Peptide / (modifications + ambig key)
    // Combine quantitated peptides - seperating out into ambiguous / unambiguous.
    Object.keys(grouped_peptides).forEach(function(pep_key) {

        // HCD identification is not used for quantification

        var peps = grouped_peptides[pep_key].filter(not_hcd_filter);
        var quan_peps = peps.filter(function(pep) { return "QuanResultID" in pep;  });
        var singlet_peps = quan_peps.filter(function(pep) { return pep.QuanChannelID.length < 2 && pep.has_pair == false; });
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
            }).map(function(pep) { return pep.areas[1] / pep.areas[0]; });

            target_ratio = Math.median( all_ratios );

            // We wish to record the mad to see if we are picking up
            // different populations of sites. Higher mad means that
            // there that there can be another population of values
            // here. This is especially useful for the ambiguous
            // modification sites, where we might see various populations.

            // We can always go back to them to examine them later.

            target_ratio_mad = Math.median(all_ratios.map(function(ratio) { return Math.abs(ratio - target_ratio); }));


            singlet_peps.forEach(function(pep) { pep.used = false; });
            ratioed_peps.forEach(function(pep) { pep.used = true; });

        } else if (singlet_peps.length > 0)  {

            target_ratio = singlet_peps[0].QuanChannelID[0] == 1 ? 1/100000 : 100000;

            var channel_ids = singlet_peps.map(function(pep) { return pep.QuanChannelID[0]; }).filter(onlyUnique);

            if (channel_ids.length > 1) {
                target_ratio = 'conflicting_singlets';
            }
        }
        if (target_ratio) {
            quan_peps.forEach(function(pep) {
                if (target_ratio_mad !== null) {
                    pep.CalculatedRatio_mad = target_ratio_mad;
                }
                pep.CalculatedRatio = target_ratio;
            });
        }
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
        var has_possible_mods = false;
        var hexnac_type = {};
        var max_score = null;

        var spectra = [];

        peps.forEach(function(pep) {
            if ("CalculatedRatio" in pep) {
                quant = pep.CalculatedRatio;
                high_sn = true;
            }
            if ("has_low_sn" in pep && ! pep.has_low_sn) {
                high_sn = true;
            }
            if ("has_pair" in pep && pep.has_pair === true && ( pep.activation !== 'HCD' && pep.activation !== 'CID' )) {
                quant = pep.QuanChannelID[0] == 1 ? 'potential_light' : 'potential_medium';
            }
            if ("hexnac_type" in pep) {
                hexnac_type[pep.hexnac_type] = true;
                if (("GlcNAc" in hexnac_type || "GalNAc" in hexnac_type)) {
                    delete hexnac_type['Unknown'];
                }
            }
            if ("possible_mods" in pep) {
                has_possible_mods = true;
            }
            spectra.push( {'score' : pep.score, 'rt' : pep.retentionTime, 'scan' : pep.scan } );
        });
        var block = {
            'multi_protein' : false,
            'sequence' : first_pep.Sequence
        };
        if (first_pep.uniprot.length > 1) {
            block.multi_protein = true;
        }
        if (quant !== null && high_sn) {
            block.quant = isNaN(parseInt(quant)) ? {'quant' : quant } : {'quant' : quant, 'mad' : first_pep.CalculatedRatio_mad };
        }
        if (Object.keys(hexnac_type).length > 0) {
            block.hexnac_type = Object.keys(hexnac_type);
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

        if (has_possible_mods) {
            block.ambiguous_mods = peps.filter(function(pep) { return pep.possible_mods; }).map(function(pep) { return write_possible_mods(pep.possible_mods); }).filter(onlyUnique);
        }
        first_pep.uniprot.forEach(function(uniprot) {
            if ( ! data[uniprot] ) {
                data[uniprot] = [];
            }
            data[uniprot].push(block);
        });
    });
    return data;
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
            return ""+mod[3]+"-"+mod[4]+"("+mod[1]+")";
        }
        return mod[0]+"("+mod[1]+")";
    }).join(',');
};


exports.combine = combine_all_peptides;