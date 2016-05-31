
var nconf = require('nconf');
var xlsx = require('node-xlsx');
var transform = require('jsonpath-object-transform');
var cellosaurus = require('./uniprot').cellosaurus;

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

var populate_source_info = function(self,organism,tissue,cell_line) {
  var results = {};
  var cell_meta;
  if (cell_line && (cell_meta = cellosaurus.lookup(cell_line))) {
    if (! cell_meta) {
      results['source-cell_line'] = cell_line;
    } else {
      results['source-cell_line'] = cell_meta.names[0];
      results['source-organism'] = cell_meta.taxid;
      results['source-organism_part'] = cell_meta.bto;
      results['source-cellosaurus_id'] = cell_meta.acc;
    }
  } else {
    results['source-organism'] = organism;
    results['source-organism_part'] = tissue;
    if (cell_line) {
      results['source-cell_line'] = cell_line;
    }
  }
  self.source = results;
  return "'source'";
};

var populate_conf = function populate_conf(manifest) {
  var manifest_data = xlsx.parse(manifest);
  if (manifest_data.length > 0) {
    var conf_data = {};
    var current_section = null;
    var current_headers = [];
    manifest_data[0].data.forEach(function(row) {
      if (row[0].match(/^\[.*\]$/)) {
        current_section = row[0].replace(/[\]\[]/g,'');
        conf_data[current_section] = conf_data[current_section] || {};
        return;
      }
      if ( ! current_section ) {
        return;
      }
      if (row[0].match(/^#/)) {
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
  var template = require('../resources/manifest_conf_mapping.template.json');
  return transform(conf_data, template, { paste: paste, summarise_ppms: summarise_ppms, populate_source_info : populate_source_info })
};


exports.populateConf = populate_conf;