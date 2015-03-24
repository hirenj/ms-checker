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

var db = new sqlite3.Database(files_to_open[0],sqlite3.OPEN_READONLY,function(err) {
	console.log(arguments);
	console.log(db);
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