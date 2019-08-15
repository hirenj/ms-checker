'use strict';

const util = require('util');
const get_cached_file = require('./ftp').get_cached_file;
const fs = require('fs');
const path = require('path');
const zlib = require('zlib');
const MultiStream = require('multistream');

// Gene Info from NCBI for selected organisms
const GENE_URLS = { fungi: 'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Fungi/All_Fungi.gene_info.gz',
                    viruses: 'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Viruses/All_Viruses.gene_info.gz',
                    inverterbrates: 'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Invertebrates/All_Invertebrates.gene_info.gz',
                    mammalia: 'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/All_Mammalia.gene_info.gz',
                    other_vertebrates : 'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Non-mammalian_vertebrates/All_Non-mammalian_vertebrates.gene_info.gz'
                  };


const get_stream = function(filename) {
  let gunzip = zlib.createGunzip();
  return fs.createReadStream(filename).pipe(gunzip);
};

// Human, CHO, Mouse, Rat, Saccharomyces Cerevisae, Saccharomyces Pombe, Baker's yeast, Drosophila, C elegans, Pig, Bovine
const valid_taxids = ['9606','10029','10090','10116','559292','4896','4932','7227','6239','9823','9913'];

const parse_files = function(stream) {
  var parse = require('csv-parse');
  var transform = require('stream-transform');

  var output = [];
  var parser = parse({delimiter: '\t', comment: '#', relax_column_count: true, quote: ''})
  var transformer = transform(function(record, callback){
    let data = record.slice(0,3);
    let taxid = data[0];
    if (valid_taxids.indexOf(taxid) < 0) {
      callback();
      return;
    }
    let geneid = data[1];
    let genename = data[2];
    let key = ['entrez',taxid,genename.toLowerCase()].join('-');
    setTimeout(function(){
      callback(null, [key,taxid,geneid,genename].join('\t')+'\n');
    }, 0);
  }, {parallel: 100});
  var stream_pipe = stream.pipe(parser).pipe(transformer).pipe(fs.createWriteStream(path.join(__dirname,'../entrezids.txt')));
  return new Promise(function(resolve,reject) {
    stream_pipe.on('end',resolve);
  });
};

const read_entries = function(stream) {
  let read_promise = Promise.resolve(true);
  try {
    fs.realpathSync(path.join(__dirname,'../entrezids.txt'));
  } catch (err) {
    read_promise = parse_files(stream);
  }
  return read_promise.then(function() {
    let reader = fs.createReadStream(path.join(__dirname,'../entrezids.txt'));
    let lineReader = require('readline').createInterface({
      input: reader
    });
    return new Promise(function(resolve,reject) {
      lineReader.on('line', function (line) {
        let data = line.split(/\t/).slice(0,4);
        let key = data[0];
        let taxid = data[1];
        let geneid = data[2];
        let genename = data[3];
        metadata[key] = (metadata[key] || []).concat( [ { entrez: geneid, name: genename, taxid: taxid } ] );
      });
      lineReader.on('close', function (line) {
        resolve();
      });
    });
  });
};

var EntrezMeta = function EntrezMeta() {

};

util.inherits(EntrezMeta,require('events').EventEmitter);

var metadata;

EntrezMeta.prototype.init = function() {
    if (metadata) {
        return Promise.resolve(true);
    }
    console.log("Initing Entrez data into ",path.join(__dirname,'../entrezids.txt'));
    let promises = Object.keys(GENE_URLS).map(function(organism) {
      return get_cached_file( GENE_URLS[organism],"entrez."+organism+".gene_info.gz", 1000*60*60*24*28 );
    });
    return Promise.all(promises).then(function(filenames) {
      metadata = {};
      if (filenames.map(function(file) { return file[1]; }).indexOf(true) >= 0) {
        try {
          fs.unlinkSync(__dirname+'/../entrezids.txt');
        } catch (err) {
          if (err.errno !== -2) {
            throw err;
          }
        }
      }
      return MultiStream(filenames.map(function(file) { return file[0]; }).map(get_stream));
    }).then(read_entries);
};

EntrezMeta.prototype.lookup = function(tax,gene) {
  return metadata[ ['entrez',tax, gene.toLowerCase() ].join('-') ];
};

exports = module.exports = new EntrezMeta();