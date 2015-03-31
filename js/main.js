// Load native UI library
var gui = require('nw.gui');

var fs 		= require('fs'),
	nconf	= require('nconf'),
	sqlite3 = require('sqlite3'),
	yauzl	= require('yauzl'),
	sax	= require('sax'),
	saxpath = require('saxpath');




var current_files = [];
nconf.env({ separator: "__", whitelist : ['MS_DEBUG'] }).overrides( require('optimist')(gui.App.argv).argv );

if (nconf.get('MS_DEBUG') || nconf.get('debug')) {
	gui.Window.get().showDevTools();
}

var files_to_open = nconf.get('_') || [];
if (files_to_open.length < 1) {
	document.getElementById('fileDialog').click();
}

var promisify_sqlite = function(db) {
	var old_all = db.all;
	db.all = function(sql,vals) {
		var args = Array.prototype.splice.call(arguments);
		return new Promise(function(resolve,reject) {
			old_all.call(db,sql,vals,function(err,vals) {
				if (err) {
					reject(err);
				} else {
					resolve(vals);
				}
			});
		});
	};
	var cached_statements = {};
	db.do_statement = function(sql,vals) {
		if ( ! cached_statements[sql]) {
			cached_statements[sql] = db.prepare(sql);
		}
		return new Promise(function(resolve,reject) {
			cached_statements[sql].all.apply( cached_statements[sql], vals.concat( function(err,vals) {
				if (err) {
					reject(err);
				} else {
					resolve(vals);
				}
			}));
		});
	};
	db.end_statement = function(sql) {
		cached_statements[sql].finalize();
		delete cached_statements[sql];
	};
};

const search_peptides_sql = 'SELECT \
	peptides.PeptideID, \
	peptides.UniquePeptideSequenceID, \
	peptides.SpectrumID, \
	peptides.Sequence, \
	PrecursorIonQuanResultsSearchSpectra.QuanResultID, \
	PrecursorIonQuanResults.QuanChannelID, \
	PrecursorIonQuanResults.Area \
FROM peptides \
	LEFT JOIN PrecursorIonQuanResultsSearchSpectra \
		ON peptides.SpectrumID = PrecursorIonQuanResultsSearchSpectra.SearchSpectrumID \
	LEFT JOIN PrecursorIonQuanResults \
		ON PrecursorIonQuanResultsSearchSpectra.QuanResultID = PrecursorIonQuanResults.QuanResultID \
WHERE peptides.ConfidenceLevel = 3 \
AND peptides.PeptideID in (SELECT distinct PeptideID \
	FROM PeptidesAminoAcidModifications \
	WHERE AminoAcidModificationID in (SELECT AminoAcidModificationID \
		FROM AminoAcidModifications \
		WHERE ModificationName like "%Hex%") \
	)';

/*
	We also want to get the subtracted HCD data. We can look at each of the FileNames, and then
	extract out all the MassPeaks (and so Spectra and Peptides) that are associated with each of
	the different filenames, and then assign the ambiguous mass difference to each of the
	peptides.
*/

const related_quants_sql = 'SELECT \
	EventID, \
	FileID, \
	Mass, \
	RT, \
	LeftRT, \
	RightRT, \
	Charge \
FROM EventAnnotations \
	JOIN Events USING(EventID) \
WHERE QuanResultID = ?';

const search_pair_sql = 'SELECT \
	* \
FROM Events \
WHERE 	FileID = ? AND \
		RT >= ? AND \
		RT <= ? AND \
		Mass >= ? AND \
		Mass <= ? \
';

const dimethyl_counts_sql = 'SELECT \
	count(distinct Position) as count, PeptideID \
FROM PeptidesAminoAcidModifications \
	JOIN AminoAcidModifications USING(AminoAcidModificationID) \
WHERE (ModificationName = "Dimethyl" OR ModificationName = "Dimethyl:2H(4)") GROUP BY PeptideID';


const peptide_metadata_sql = 'SELECT \
	Description \
FROM Peptides \
	LEFT JOIN PeptidesProteins USING(PeptideID) \
	LEFT JOIN ProteinAnnotations USING (ProteinID) \
WHERE Peptides.PeptideID = ?';

// const peptide_modification_sql = 'SELECT \
// 	Position, ModificationName \
// FROM PeptidesAminoAcidModifications \
// 	LEFT JOIN AminoAcidModifications USING (AminoAcidModificationID) \
// WHERE PeptideID = ?';

const peptide_modification_sql = 'SELECT \
	Position,AminoAcidModificationID \
FROM PeptidesAminoAcidModifications \
WHERE PeptideID = ?';

const all_peptide_modifications_sql = 'SELECT \
	PeptideID, Position, ModificationName \
FROM PeptidesAminoAcidModifications \
	LEFT JOIN AminoAcidModifications USING (AminoAcidModificationID)';

const retrieve_spectrum_sql = 'SELECT \
	Spectrum \
FROM SpectrumHeaders \
	LEFT JOIN Spectra USING (UniqueSpectrumID) \
WHERE SpectrumID = ?';

const MASS_MEDIUM = 32.056407;
const MASS_LIGHT = 28.0313;

var dimethyl_count_cache = null;

var find_dimethyls = function(db,pep) {
	return new Promise(function(resolve,reject) {

		if (! dimethyl_count_cache ) {

			dimethyl_count_cache = {};
			console.log("Populating Dimethyl cache");

			db.all(dimethyl_counts_sql,[]).then(function(data) {
				dimethyl_count_cache = {};
				data.forEach(function(count) {
					dimethyl_count_cache[count.PeptideID] = count.count+1;
				});
				console.log("Populated Dimethyl cache");

				if ( ! pep.QuanResultID ) {
					resolve();
				} else {
					resolve(find_dimethyls(db,pep));
				}
			},reject);

		} else {
			resolve(dimethyl_count_cache[pep.PeptideID]);
			// return check_potential_pair(db,pep,dimethyl_count_cache[pep.PeptideID]);
		}
	});
};


var validated_quans_cache = {};

var check_potential_pair = function(db,pep,num_dimethyl) {
	return new Promise(function(resolve,reject) {
		if (! pep.QuanResultID) {
			resolve();
			return;
		}
		if (pep.QuanResultID in validated_quans_cache) {
			pep.has_pair = validated_quans_cache[pep.QuanResultID];
			resolve();
			return;
		}
		if (pep.Sequence.indexOf('K') < 0 ) {
			num_dimethyl = 1;
		}
		if ( ! num_dimethyl ) {
			find_dimethyls(db,pep).then(function(count) {
				resolve(check_potential_pair(db,pep,count));
			},reject);
			return;
		}

		pep.has_pair = false;

		var mass_change_dir = 1;

		if (pep.QuanChannelID.indexOf(2) >= 0) {
			mass_change_dir = -1;
		}

		db.all(related_quants_sql,[ pep.QuanResultID ]).then(function(events) {
			var events_length = events.length;
			events.forEach(function(ev) {
				events_length -= 1;
				if (pep.has_pair) {
					if (events_length <= 0) {
						validated_quans_cache[pep.QuanResultID] = pep.has_pair;
						resolve();
					}
					return;
				}
				var target_mass = ev.Mass + (mass_change_dir * num_dimethyl * (MASS_MEDIUM-MASS_LIGHT)/ev.Charge);
				db.all(search_pair_sql, [ ev.FileID, ev.LeftRT - 0.05, ev.RightRT + 0.05, target_mass*(1 - 15/1000000) , target_mass*(1 + 15/1000000)]).then(function(rows) {
					if (rows && rows.length > 0) {
						pep.has_pair = true;
						if (nconf.get('MS_DEBUG') || nconf.get('debug')) {
							pep.matching_events = rows.map(function(row) { return row.EventID; });
							pep.search_event = ev.EventID;
							pep.target_mass = target_mass;
							pep.num_dimethyl = num_dimethyl;
							pep.mass_change_dir = mass_change_dir;
						}
					}
					if (events_length <= 0) {
						validated_quans_cache[pep.QuanResultID] = pep.has_pair;
						resolve();
					}
				}).catch(reject);
			});
		}).catch(reject);
	});
};

var pep_calls = 0;

var produce_peptide_data = function(db,pep) {
	return db.do_statement(peptide_metadata_sql, [ pep.PeptideID ]).then(function(pep_datas) {
		pep_calls++;
		if ((pep_calls % 100) == 0) {
			console.log(pep_calls);
		}
		if ( ! pep.uniprot ) {
			pep.uniprot = [];
		}
		pep_datas.forEach(function(pep_data) {
			var uniprot = pep_data.Description.split('|')[1];
			if (pep.uniprot.indexOf(uniprot) < 0) {
				pep.uniprot.push(uniprot);
			}
		});
	});
};

var peptide_modifications_cache = null;

var produce_peptide_modification_data = function(db,pep) {
	return new Promise(function(resolve,reject) {

		if (! peptide_modifications_cache ) {

			peptide_modifications_cache = {};
			console.log("Populating Modifications cache");

			db.all(all_peptide_modifications_sql,[]).then(function(data) {
				peptide_modifications_cache = {};
				data.forEach(function(mods) {
					peptide_modifications_cache[mods.PeptideID] =  peptide_modifications_cache[mods.PeptideID] || [];
					peptide_modifications_cache[mods.PeptideID].push([ mods.Position + 1, mods.ModificationName ]);
				});
				console.log("Populated Modifications cache");

				if ( ! pep.PeptideID ) {
					resolve();
				} else {
					resolve(produce_peptide_modification_data(db,pep));
				}
			},reject);

		} else {
			pep.modifications = (peptide_modifications_cache[pep.PeptideID] || []).filter( function(mod) { return (mod[1] || '').indexOf('Hex') >= 0 } ).map( function(mod) { return mod } );
			resolve();
		}
	});
};

var combine_peptides = function(peps) {
	var all_peps = {};
	peps.forEach(function(pep) {
		var curr_pep = all_peps[pep.PeptideID] ||  pep;
		curr_pep.areas = curr_pep.areas || [];

		if (! Array.isArray(curr_pep.QuanChannelID) ) {
			curr_pep.QuanChannelID =  curr_pep.QuanChannelID ? [ curr_pep.QuanChannelID ] : [];
		}
		if (curr_pep != pep && pep.QuanChannelID) {
			curr_pep.QuanChannelID.push(pep.QuanChannelID);
		}
		if (pep.Area) {
			curr_pep.areas.push(pep.Area);
		}
		all_peps[pep.PeptideID] = curr_pep;
	});
	return Object.keys(all_peps).map(function(key) { return all_peps[key]; });
};

var check_spectrum = function(db,pep) {

	var saxParser = sax.createStream(true)
	var streamer = new saxpath.SaXPath(saxParser, '//PeakCentroids/Peak');
	var streamer2 = new saxpath.SaXPath(saxParser, '//ScanEvent/ActivationTypes');

	streamer.on('match', function(xml) {
		// Check to see what the ratio (mz-138 + mz-168) / (mz-126 + mz-144) is
		// if it is within 0.4 - 0.6, it is a GalNAc, 2.0 or greater, GlcNAc
	    console.log(xml);
	});

	streamer2.on('match', function(xml) {
		// Check that we're matching a HCD spectrum here
	    console.log(xml);
	});

	return db.all(retrieve_spectrum_sql, [ pep.SpectrumID ]).then(function(spectra) {
		var spectrum = spectra[0].Spectrum;

		yauzl.fromBuffer(spectrum, function(err, zipfile) {

		  if (err) throw err;

		  zipfile.on("entry", function(entry) {
		    if (/\/$/.test(entry.fileName)) {
		      // directory file names end with '/'
		      return;
		    }
		    zipfile.openReadStream(entry, function(err, readStream) {
		      if (err) throw err;
		      readStream.pipe(saxParser);
		    });
		  });
		});
	});
};


var onlyUnique = function(value, index, self) {
    return self.indexOf(value) === index;
};

var global_results;

var db = new sqlite3.Database(files_to_open[0],sqlite3.OPEN_READONLY,function(err) {

	promisify_sqlite(db);
	produce_peptide_modification_data(db,{}).then(function() { return find_dimethyls(db,{'PeptideID': 0}); }).then(function() {
		console.log("Searching Peptides");
		db.all(search_peptides_sql).then(function(peps) {
			console.log("Retrieved peptides to search: "+(peps.length)+" total peptides");

			var combined_peps = combine_peptides(peps);

			var singlet_peps = combined_peps.filter( function(pep) { return pep.QuanChannelID.filter(onlyUnique).length == 1; } );

			var pair_promises = singlet_peps.map( function(pep) { return check_potential_pair(db,pep,null); } );

			var metadata_promises = combined_peps.map( function(pep) { return produce_peptide_data(db,pep).then( produce_peptide_modification_data(db,pep) );  } );

			Promise.all(pair_promises.concat(metadata_promises)).then(function() {
				db.end_statement(peptide_metadata_sql);
				console.log("Done");
				global_results = combined_peps; //combined_peps.filter(function(pep) { return pep.QuanChannelID.filter(onlyUnique).length == 1; }) );
				global_db = db;
			}, function() { console.log("Rejected"); });
		});
	});


	// Area is the sum of Areas for the Events associated with the QuanResultID (from EventAnnotations)

	// select PeptideId,SpectrumID from peptides where peptides.ConfidenceLevel = 3

	// select QuanResultID, SearchSpectrumID from PrecursorIonQuanResultsSearchSpectra where SearchSpectrumID = ?

	// select * from PrecursorIonQuanResults where QuanResultID = ?

	// We want to search:
	//  This looks like this is all the spectra that were used to search for an area
	//	PrecursorIonAreaSearchSpectra


	//  This is just the search spectra where an area was able to be obtained.
	//	PrecursorIonQuanResults	
	//	PrecursorIonQuanResultsSearchSpectra
	//  We want to look at the SearchSpectrumID to get the Peptide / Spectrum link

	// CREATE TABLE [PrecursorIonQuanResults] (
	// 	[QuanChannelID][int] NOT NULL,
	// 	[QuanResultID][int] NOT NULL,
	// 	[Mass] [double] NOT NULL ,
	// 	[Charge] [int] NOT NULL ,
	// 	[Area] [double] NULL ,
	// 	[RetentionTime] [double] NULL
	// );
	// CREATE TABLE [PrecursorIonQuanResultsSearchSpectra] (
	// 	[ProcessingNodeNumber] [int] NOT NULL ,
	// 	[QuanResultID] [int] NOT NULL ,
	// 	[SearchSpectrumID] [int] NULL 
	// );

	// CREATE TABLE [Peptides] (	
	// 	[ProcessingNodeNumber] [int] NOT NULL ,
	// 	[PeptideID] [int] NOT NULL ,
	// 	[SpectrumID] [int] NOT NULL ,
	// 	[TotalIonsCount] [smallint] NOT NULL ,
	// 	[MatchedIonsCount] [smallint] NOT NULL ,
	// 	[ConfidenceLevel] [smallint] NOT NULL ,
	// 	[SearchEngineRank] [int] NOT NULL,
	// 	[Hidden] [bit] DEFAULT 0 NOT NULL ,
	// 	[Sequence] [varchar] (1024) COLLATE NOCASE NULL,	
	// 	[Annotation] [varchar] (1024) COLLATE NOCASE NULL ,
	// 	[UniquePeptideSequenceID] [int] DEFAULT 1 NOT NULL,
	// 	[MissedCleavages] [smallint] NOT NULL,
	// 	PRIMARY KEY (		
	// 		[ProcessingNodeNumber],
	// 		[PeptideID] 
	// 	)
	// );

	// EventAnnotations and EventAreaAnnotations should provide the link
	// between the EventID and QuanResultID for a peptide

	// CREATE TABLE [EventAnnotations] (
	// 	[EventID] INTEGER PRIMARY KEY NOT NULL,
	// 	[Charge] [smallint] NOT NULL,									 
	// 	[IsotopePatternID] [int] NOT NULL,
	// 	[QuanResultID] [int] NOT NULL,
	// 	[QuanChannelID] [int] NOT NULL
	// );

	// CREATE TABLE [EventAreaAnnotations] (
	// 	[EventID] INTEGER PRIMARY KEY NOT NULL,
	// 	[Charge] [smallint] NOT NULL,									 
	// 	[IsotopePatternID] [int] NOT NULL,
	// 	[QuanResultID] [int] NOT NULL
	// );

});