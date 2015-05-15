var xlsx = require('node-xlsx');
var fs = require('fs');

var write_excel_file = function(datablock,filename) {
	var rows = [[ "uniprot", "peptide_id", "quant" , "quant_mad", "hexnac_type", "sequence", "composition", "spectra",  "site", "ambiguous"  ]];
	var peptide_id = 0;
    Object.keys(datablock).forEach(function(uniprot) {
    	var peps = datablock[uniprot].filter(function(pep) { return ! pep.multi_protein; });
    	peps.forEach(function(pep) {
	    	var data = [uniprot];

    		peptide_id += 1;
    		data.push(peptide_id);

    		if (pep.quant) {
    			data.push(pep.quant.quant);
    			data.push(pep.quant.mad);
    		} else {
    			data.push(null);
    			data.push(null);
    		}
    		if (pep.hexnac_type) {
    			data.push(pep.hexnac_type);
    		} else {
    			data.push(null);
    		}
    		data.push(pep.sequence);
            data.push(pep.composition.join(';'));
            data.push(pep.spectra.map(function(spec) {  return spec.score +"-" + spec.rt + "-" + spec.scan;  }).join(","));

    		if (pep.sites) {
    			pep.sites.forEach(function(site) {
    				rows.push(data.concat([site[0]]));
    			});
    		}
    		if (pep.ambiguous_mods) {
    			pep.ambiguous_mods.forEach(function(ambig) {
    				rows.push(data.concat([null,ambig]));
    			});    			
    		}
    	});
    });
    return new Promise(function(resolve,reject) {
        fs.writeFile( filename, xlsx.build([{name:"Peptides", data: rows }]) , function(err) {
            if (err) {
                reject(err);
            }
            resolve();
        });
    });
};

exports.write = write_excel_file;