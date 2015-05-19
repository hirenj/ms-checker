var xlsx = require('node-xlsx');
var fs = require('fs');

var write_excel_file = function(datablock,filename) {
    var rows = [[ "source", "uniprot", "peptide_id", "quant" , "quant_mad", "hexnac_type", "sequence", "composition", "spectra",  "site", "ambiguous"  ]];
    var metadata = [];
    var peptide_id = 0;
    if (! Array.isArray(datablock.metadata) && datablock.metadata ) {
        datablock.metadata = [ datablock.metadata ];
    }
    (datablock.metadata || []).forEach(function(meta) {
        metadata.push([ JSON.stringify(meta,null,'    ') ]);
    });
    Object.keys(datablock.data).forEach(function(uniprot) {
        var peps = datablock.data[uniprot].filter(function(pep) { return ! pep.multi_protein; });
        peps.forEach(function(pep) {

            var data = [];
            data.push(pep.source || null);
            data.push(uniprot);

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
            if (! pep.sites && ! pep.ambiguous_mods) {
                rows.push(data.concat([ null,null ]));
            }
        });
    });
    return new Promise(function(resolve,reject) {
        fs.writeFile( filename, xlsx.build([{name:"Peptides", data: rows },{name: "Metadata", data: metadata }]) , function(err) {
            if (err) {
                reject(err);
            }
            resolve();
        });
    });
};

exports.write = write_excel_file;