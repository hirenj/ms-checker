'use strict';

var nconf = require('nconf');
var xlsx = require('node-xlsx');
var transform = require('jsonpath-object-transform');
var cellosaurus = require('./uniprot').cellosaurus;
var entrez = require('./entrez');
var Validator = require('jsonschema').Validator;

var paste = function(self,sep) {
  var args = Array.prototype.slice.call(arguments);
  var cols = args.slice(2);
  var arrays = cols.map(function(col) { return self[col] });
  self.pasted = arrays.reduce(function(curr,next) {
    if ( ! curr ) {
      return next;
    }
    return curr.map(function(val,idx) { return val + sep + next[idx]; });
  });
  return "'pasted'";
};

var summarise_ppms = function(self,ppm,sources,filenames) {
  var results = {};
  filenames.map(function(filename,i) {
    var pattern_match = ppm.MSF_pattern.map(function(pattern) { return pattern === '*' || filename.indexOf(pattern) >= 0;  });
    var ppm_index = pattern_match.indexOf(true);
    if (ppm_index < 0) {
      return;
    }
    results[ "ppm-"+sources[i]+"-min" ] = ppm.Minimum[ppm_index];
    results[ "ppm-"+sources[i]+"-max" ] = ppm.Maximum[ppm_index];
  });
  self.parsed_ppm = results;
  return "'parsed_ppm'";
};

var summarise_quants = function(self,channels,identifiers) {
  var results = {};
  channels.forEach(function(channel,i) {
    results[ "quant-channel-"+channel ] = identifiers[i];
  });
  self.parsed_quants = results;
  return "'parsed_quants'";
};


var summarise_flags = function(self) {
  var results = {};
  if (self.Flags) {
    self.Flags.Flag.forEach(function(flag,i) {
      results[ flag ] = self.Flags.Value[i];
    });
  }
  self.flags = results;
  return "'flags'";
};

var populate_source_info = function(self,organism,tissue,cell_line) {
  var results = {};
  var cell_meta;
  if (cell_line && (cell_meta = cellosaurus.lookup(cell_line))) {
    if (! cell_meta) {
      results['source-cell_line'] = cell_line;
    } else {
      results['source-cell_line'] = cell_meta.names[0];
      results['source-organism'] = cell_meta.taxid;
      results['source-organism_part'] = cell_meta.bto ? 'bto:'+cell_meta.bto : tissue[0];
      results['source-cellosaurus_id'] = cell_meta.acc;
    }
  } else {
    results['source-organism'] = organism[0];
    results['source-organism_part'] = tissue[0];
    if (cell_line) {
      results['source-cell_line'] = cell_line;
    }
  }
  self.source = results;
  return "'source'";
};

var populate_perturbation_info = function(self,taxid,genes,types,sample_identifier) {

  sample_identifier = sample_identifier.map( ids => ids.split(',').map( id => id.trim().toLowerCase() ) );

  let new_result = () => {return { 'perturbation-ko': [], 'perturbation-ki': [], 'perturbation-wt' : [] } };
  let all_samples = Array.concat.apply([],sample_identifier).filter( id => id !== 'wt');
  let all_results = {};
  for (let sample of ['wt'].concat(all_samples)) {
    let result = new_result();
    all_results[sample] = result;
    genes.forEach(function(gene,idx) {
      // Skip genes that do not have the sample identifier
      // If there's no identifier, we assume it is the common sample (wt)
      if ( sample_identifier[idx].indexOf(sample) < 0 && sample_identifier[idx].length > 0 ) {
        return;
      }
      let gene_data = entrez.lookup(taxid,gene)[0];
      let type = 'perturbation-'+types[idx].toLowerCase();
      if ( ! result[type] ) {
        result[type] = [];
      }
      result[type].push({'entrez' : parseInt(gene_data.entrez), 'symbol' : gene_data.name })
    });
    if (self.Other_Perturbations && self.Other_Perturbations.Description.length > 0 && (! self.Other_Perturbations.Identifier || self.Other_Perturbations.Identifier.toLowerCase().indexOf(sample) >= 0)) {
      result['perturbation-other'] = self.Other_Perturbations.Description[0];
    } else {
      result['perturbation-other'] = [];
    }
  }
  self.perturbation = all_samples.length === 1 ? all_results[all_samples[0]] : all_results;
  return "'perturbation'";
};

var populate_conf = function populate_conf(manifest) {
  var manifest_data = xlsx.parse(manifest);
  if (manifest_data.length > 0) {
    var conf_data = {};
    var current_section = null;
    var current_headers = [];
    manifest_data[0].data.forEach(function(row) {
      if (row[0] && (row[0]+'').match(/^\[.*\]$/)) {
        current_section = row[0].replace(/[\]\[]/g,'').replace(/\s/g,'_');
        conf_data[current_section] = conf_data[current_section] || {};
        return;
      }
      if ( ! current_section ) {
        return;
      }
      if (row[0] && (row[0]+'').match(/^#/)) {
        current_headers = row.map(function(header) { return header.replace('#','').replace(/\s/g,'_'); });
        current_headers.forEach(function(header) {
          conf_data[current_section][header] = [];
        });
        return;
      }
      row.forEach(function(data,i) {
        conf_data[current_section][current_headers[i]].push(row[i]);
      });
    });
  }
  var manifest_version = conf_data.Manifest.Version[0];
  if ( ! manifest_version.match(/^\d+\.?\d+?$/)) {
    throw new Error('Invalid manifest');
  }
  var schema = require('../resources/manifests/'+manifest_version+'/schema.json');
  var template = require('../resources/manifests/'+manifest_version+'/manifest_conf_mapping.template.json');
  var valid = (new Validator()).validate(conf_data,schema);
  if ( valid.errors && valid.errors.length > 0 ) {
    throw new Error(valid.errors);
  }
  console.log("Validated manifest");
  return transform(conf_data, template, { paste: paste,
                                          summarise_ppms: summarise_ppms,
                                          summarise_quants : summarise_quants,
                                          summarise_flags : summarise_flags,
                                          populate_source_info: populate_source_info,
                                          populate_perturbation_info: populate_perturbation_info
                                        })
};


exports.populateConf = populate_conf;