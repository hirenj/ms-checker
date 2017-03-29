var fs      = require('fs'),
    sqlite3 = require('sqlite3');

var nconf = require('nconf');
var ProgressBar = require('progress');

Math.median = require('../js/math-median').median;

var uniprot_meta = require('../js/uniprot.js').uniprot;

var promisify_sqlite = require('../js/sqlite-promise');

var quantitative = require('../js/quantitative');
var ambiguous = require('../js/ambiguous');
var hexnac_hcd = require('../js/hexnac_hcd');
var fragmentation = require('../js/fragmentions');
var metadata = require('../js/metadata');
var ppm = require('../js/ppm');

var peptide = require('../js/peptide');

var spectra = require('../js/spectrum');

var contaminants = require('../js/contaminants');

var ProgressBar = require('progress');

ambiguous.conf = nconf;
hexnac_hcd.conf = nconf;

var progress_bars = {};

function Multibar(stream) {
  this.stream     = stream || process.stderr;
  this.cursor     = 0;
  this.bars       = [];
  this.terminates = 0;
};

Multibar.prototype = {
  newBar: function(schema, options) {
    options.stream = this.stream;
    var bar = new ProgressBar(schema, options);
    this.bars.push(bar);
    var index = this.bars.length - 1;

    // alloc line
    this.move(index);
    this.stream.write('\n');
    this.cursor ++;

    // replace original
    var self  = this;
    bar.otick = bar.tick;
    bar.oterminate = bar.terminate;
    bar.tick = function(value, options) {
      self.tick(index, value, options);
    }
    bar.terminate = function() {
      self.terminates++;
      if (self.terminates == self.bars.length) {
        self.terminate();
      }
    }

    return bar;
  },

  terminate: function() {
    this.move(this.bars.length);
    this.stream.clearLine();
    this.stream.cursorTo(0);
  },

  move: function(index) {
    if (!this.stream.isTTY) return;
    this.stream.moveCursor(0, index - this.cursor);
    this.cursor = index;
  },

  tick: function(index, value, options) {
    var bar = this.bars[index];
    if (bar) {
      this.move(index);
      bar.otick(value, options);
    }
  }
};

var mbars = new Multibar();

[quantitative,ambiguous,hexnac_hcd,fragmentation,peptide,ppm].forEach(function(module) {
    module.on('task',function(desc) {
        if ( !process.stdin.isTTY ) {
          var max = -1;
          this.addListener('progress',function(percentage) {
              var barname = (module.constructor.name+" "+desc+(Array(50).join(' '))).substring(0,50);
              var percentage_level = parseInt((10 * percentage).toFixed(0));
              if (percentage_level > max) {
                  max = percentage_level;
                  console.log(barname + " " + max + "0 %");
                  if (percentage_level == 10) {
                    module.removeListener('progress',arguments.callee);
                  }
              }
          });
          return;
        }
        this.addListener('progress',function(percentage) {
            if ( ! progress_bars[module.constructor.name+desc] ) {
                var barname = (module.constructor.name+" "+desc+(Array(50).join(' '))).substring(0,50);
                progress_bars[module.constructor.name+desc] = mbars.newBar(barname+' [:bar] ', { total: 100} );
            }

            var progress = progress_bars[module.constructor.name+desc];
            progress.tick( (100*percentage).toFixed(0) - progress.curr );
            if (percentage == 1) {
                delete progress_bars[module.constructor.name+desc];
                module.removeListener('progress',arguments.callee);
            }
        });
    });
});


var partial = function(fn) {
    var args = Array.prototype.slice.call(arguments, 1);
    return function() { return fn.apply(this, args.concat(Array.prototype.slice.call(arguments, 0))); };
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

var open_db = function(filename) {
    return new Promise(function(resolve,reject) {
        var db = new sqlite3.Database(filename,sqlite3.OPEN_READWRITE,function(err) {
            if (err) {
              console.log(err);
              reject(err);
              return;
            }
            db.run('CREATE index if not exists ms_check_spectrum_peak_id on SpectrumHeaders(MassPeakID)');
            if (err) {
                reject(err);
                return;
            }
            resolve(db);
        });
    });
};

var browse_file = function(filename) {
    return open_db(filename).then(function(db) {
        promisify_sqlite(db);
        return db;
    });
};

var tracker_labels = {};

var tracker = function(label) {
  if (nconf.get('track-peptide-sequence')) {
    return (function(peps) {
      var filtered = peps.filter(function(pep) { return pep.Sequence == nconf.get('track-peptide-sequence'); });
      if (filtered.length > 0) {
//        console.log(filtered);
      } else {
        if ( ! tracker_labels[nconf.get('track-peptide-sequence')] ) {
          console.log("\n\nWe lost tracked peptide "+nconf.get('track-peptide-sequence')+" at "+label+"\n\n");
          tracker_labels[nconf.get('track-peptide-sequence')] = true;
        }
      }
      return peps;
    });
  }
  if (nconf.get('track-quanid')) {
    return (function(peps) {
      var filtered = peps.filter(function(pep) { return pep.QuanResultID == nconf.get('track-quanid'); });
      if (filtered.length > 0) {
//        console.log(filtered);
      } else {
        if ( ! tracker_labels[nconf.get('track-quanid')] ) {
          console.log("\n\nWe lost tracked peptide "+nconf.get('track-quanid')+" at "+label+"\n\n");
          tracker_labels[nconf.get('track-quanid')] = true;
        }
      }
      return peps;
    });
  }
  if (nconf.get('track-peptide-count')) {
    return (function(peps) {
      if (peps.length === 0 && ! tracker_labels['any_peptides']) {
        console.log("No peptides left at ",label);
        tracker_labels['any_peptides'] = true;
      }
      return peps;
    });
  }
  return function(peps) { return peps; }
};

var process_data = function(filename,sibling_files,source) {
    var global_datablock = {'data' : {}, 'metadata' : {}};
    return open_db(filename).then(function(db) {
        db.source = source;
        promisify_sqlite(db);
        metadata.get_metadata(db).then(function(meta) {
            global_datablock.metadata = meta;
        });
        var calculated_cutoffs = [];

        var ppms_promise = ppm.select_ppm_cutoffs(db).then(function(ppm_cutoffs) {
          ppm_cutoffs.forEach(function(cutoff) {
            calculated_cutoffs.push(cutoff);
          });
        });

        var peptides_promise = ppms_promise.then(function() {

          var quant_promise = quantitative.init_caches(db).then(function() {
              return quantitative.retrieve_quantified_peptides(db).then(tracker('Retrieve quantified peptides')).then(partial(quantitative.check_quantified_peptides,db)).then(tracker('Check quantified peptides'));
          });

          var ambig_promise = Promise.resolve(true).then(function() {
              return ambiguous.retrieve_peptides(db).then(function(peps) { return peps? peps : []; });
          });

          return Promise.all([quant_promise,ambig_promise]);
        });

        var processing_promise = peptides_promise.then(function(all_peps) {
            var merged = Array.prototype.concat.apply([], all_peps);
            return merged;
        }).then(tracker('Read from DB'))
          .then(peptide.produce_peptide_scores_and_cleanup.bind(peptide,db))
          .then(tracker('Scores'))
          .then(peptide.calculate_deltacn.bind(peptide))
          .then(peptide.filter_deltacn.bind(peptide,0.05))
          .then(tracker('Delta CN filtering'))
          .then(peptide.filter_ambiguous_spectra)
          .then(tracker('Filter ambiguous spectra'))
          .then(hexnac_hcd.guess_hexnac.bind(hexnac_hcd,db,sibling_files))
          .then(tracker('HexNAc HCD'))
          .then(quantitative.populate_missing_quant_channel)
          .then(tracker('Populate missing quant channel'))
          .then(partial(fragmentation.validate_peptide_coverage,db))
          .then(tracker('Fragmentation validation'))
          .then(peptide.merge_modifications_deltacn.bind(peptide))
          .then(tracker('Delta CN based modification merging'))
          .then(ppm.annotate_ppm)
          .then(tracker('PPM annotation'))
          .then(ppm.filter_peptides.bind(ppm,calculated_cutoffs))
          .then(spectra.filter_hcd_with_mods)
          .then(tracker('HCD modification filtering'));


        return processing_promise.then(function(peps) {
            tracker('Finished tracking, did not lose')([]);
            return uniprot_meta.init().then(function() {
                global_datablock.data = peptide.combine(peps);
                quantitative.clear_caches();
                spectra.clear_caches();
                fragmentation.clear_caches();
                return global_datablock;
            });
        });
    });
};


var files_by_source = {};

var process_files = function(files,sources) {
    var result = Promise.resolve(true);
    var blocks = [];
    files = [].concat(files);
    sources = [].concat(sources);
    files.forEach(function(file,idx) {
        var source = sources[idx];
        files_by_source[source] = (files_by_source[source] || []).concat(file);
    });
    return result.then(function() {
        var self_func = arguments.callee;
        var file = files.shift();
        var source = sources.shift();
        if ([].concat(nconf.get('option-swap-channels')).indexOf(source) >= 0) {
          quantitative.swapChannelLabels(true);
        } else {
          quantitative.swapChannelLabels(false);
        }


        if (! file) {
            return blocks;
        }
        var other_files = files_by_source[source].filter(function(other_file) {
            return other_file !== file;
        });
        return process_data(file,other_files,source).then(function(block) {
            blocks.push(block);
            return self_func();
        });
    });
};

var combine = function(blocks,sources) {
    var result = { data : {} , metadata : [] };
    blocks.forEach(function(block) {
        var source = (sources || []).shift();

        Object.keys(block.data).forEach(function(prot) {
            if ( ! result.data[prot] ) {
                result.data[prot] = [];
            }
            if (source) {
                block.data[prot] = block.data[prot].map(function(pep) {
                    return clone(pep);
                });
                block.data[prot].forEach(function(pep) {
                    pep.source = source;
                });
            }
            result.data[prot] = result.data[prot].concat(block.data[prot]);
        });
        result.metadata.push(block.metadata);
    });

    result.metadata = metadata.combine(result.metadata);
    result.metadata.title = nconf.get('output').replace('.msdata','');
    return contaminants.get_version({}).then(function() {
      contaminants.identifiers.forEach(function(contaminant) {
        if (contaminant in result.data) {
          delete result.data[contaminant];
        }
      });
      delete result.data[undefined];
      return result;
    });
};

module.exports = exports = { process: process_files, combine: combine, browse_file: browse_file };

