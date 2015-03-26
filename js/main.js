// Load native UI library
var gui = require('nw.gui');

var fs 		= require('fs'),
	nconf	= require('nconf'),
	sqlite3 = require('sqlite3');




var current_files = [];
nconf.env({ separator: "__", whitelist : ['MS_DEBUG'] }).overrides( require('optimist')(gui.App.argv).argv );
console.log(nconf.get());

if (nconf.get('MS_DEBUG') || nconf.get('debug')) {
	gui.Window.get().showDevTools();
}

var files_to_open = nconf.get('_') || [];
if (files_to_open.length < 1) {
	document.getElementById('fileDialog').click();
}

var sql_string = 'SELECT \
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

var related_quants = 'SELECT \
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

var search_pair_sql = 'SELECT \
	* \
FROM Events \
WHERE 	FileID = ? AND \
		RT >= ? AND \
		RT <= ? AND \
		Mass >= ? AND \
		Mass <= ? \
';

var dimethyl_counts = 'SELECT \
	count(distinct Position) as count, PeptideID \
FROM PeptidesAminoAcidModifications \
	JOIN AminoAcidModifications USING(AminoAcidModificationID) \
WHERE (ModificationName = "Dimethyl" OR ModificationName = "Dimethyl:2H(4)") GROUP BY PeptideID';


var dimethyl_count_cache = null;

var find_dimethyls = function(db,pep,callback) {
	if (! dimethyl_count_cache ) {
		dimethyl_count_cache = {};
		console.log("Populating Dimethyl cache");
		db.all(dimethyl_counts,[],function(err,data) {
			dimethyl_count_cache = {};
			data.forEach(function(count) {
				dimethyl_count_cache[count.PeptideID] = count.count;
			});
			setTimeout(function() {
				console.log("Populated Dimethyl cache");
				find_dimethyls(db,pep,callback);
				if ( ! pep.QuanResultID ) {
					callback();
				}
			},0);
		});
		return;
	}
	setTimeout(function() {
		check_potential_pair(db,pep,dimethyl_count_cache[pep.PeptideID],callback);
	},0);
};



var validated_quans = {};

var check_potential_pair = function(db,pep,num_dimethyl,callback) {
	if (! pep.QuanResultID) {
		return;
	}
	if (pep.QuanResultID in validated_quans) {
		pep.has_pair = validated_quans[pep.QuanResultID];
		callback();
		return;
	}
	if (pep.Sequence.indexOf('K') < 0 ) {
		num_dimethyl = 1;
	}
	if ( ! num_dimethyl ) {
		find_dimethyls(db,pep,callback);
		return;
	}
	pep.has_pair = false;
	var mass_change_dir = 1;
	if (pep.QuanChannelID.indexOf(2) >= 0) {
		mass_change_dir = -1;
	}
	db.all(related_quants,[ pep.QuanResultID ],function(err,events) {
		var events_length = events.length;
		events.forEach(function(ev) {
			events_length -= 1;
			if (pep.has_pair) {
				if (events_length <= 0) {
					validated_quans[pep.QuanResultID] = pep.has_pair;
					callback();
				}
				return;
			}
			var target_mass = ev.Mass + (mass_change_dir * num_dimethyl * (32.056407-28.0313)/ev.Charge);
			// console.log(ev.Mass);
			// console.log(ev.Charge);
			// console.log(pep.QuanChannelID);
			// console.log(target_mass);
			db.all(search_pair_sql, [ ev.FileID, ev.LeftRT - 0.05, ev.RightRT + 0.05, target_mass*(1 - 15/1000000) , target_mass*(1 + 15/1000000)],function(err,rows) {
				if (rows && rows.length > 0) {
					pep.has_pair = true;	
				}
				if (events_length <= 0) {
					validated_quans[pep.QuanResultID] = pep.has_pair;
					callback();
				}
			});
		});
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

var db = new sqlite3.Database(files_to_open[0],sqlite3.OPEN_READONLY,function(err) {
	var to_validate = [];
	var validated = [];

	find_dimethyls(db,{'PeptideID': 0},function() {
		console.log("Searching Peptides");
		db.all(sql_string,function(err,peps) {
			console.log("Retrieved peptides to search: "+(peps.length)+" total peptides");
			singlet_peps = combine_peptides(peps).filter(function(pep) {
				return pep.QuanChannelID.length == 1;
			});
			singlet_peps.forEach(function(pep) {
				check_potential_pair(db,pep,null,function(masses) {
					if (pep.has_pair) {
						if (to_validate.indexOf(pep.QuanResultID) < 0) {
							to_validate.push(pep.QuanResultID);
							console.log(pep);
						}
						console.log("Need to validate "+to_validate.length);
					} else {
						if (validated.indexOf(pep.QuanResultID) < 0) {
							validated.push(pep.QuanResultID);
						}
						console.log("Now validated "+validated.length);
					}
				});
			});
			// console.log(arguments);
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