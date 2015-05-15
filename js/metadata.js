var fs = require('fs');

const fasta_files_sql = "SELECT FileName FROM FastaFiles";

const raw_files_sql = "SELECT FileName, FileTime FROM FileInfos";

const pd_metadata_sql = "SELECT SoftwareVersion, Date FROM SchemaInfo";

const score_type_sql = 'SELECT ScoreID, ScoreName, Description FROM ProcessingNodeScores';

const MSDATA_FORMAT_VERSION = "1";

var get_raw_filenames = function(db,metadata){
    var raw_filenames = [];
    return db.each(raw_files_sql,[],function(err,file) {
        raw_filenames.push({'file': file.FileName, 'created': new Date(file.FileTime).toISOString() });
    }).then(function() {
        metadata['ms_run-location'] = raw_filenames;
    });
};

var get_fasta_filenames = function(db,metadata) {
    var fasta_filenames = [];
    return db.each(fasta_files_sql,[],function(err,file) {
        fasta_filenames.push(file.FileName);
    }).then(function() {
        metadata['database'] = fasta_filenames;
        metadata['database_version'] = fasta_filenames.map(function() {  return "" });
    });
};

var get_pd_metadata = function(db,metadata) {
    var version;
    var run_date;
    return db.each(pd_metadata_sql,[],function(err,software) {
        version = software.SoftwareVersion;
        run_date = software.Date;
    }).then(function() {
        metadata['software'] = metadata['software'] || [];
        metadata['software'].push({
            'version'   : version,
            'name'      : 'Proteome Discoverer',
            'run-date'  : run_date
        });
    });
};


var get_self_version = function(metadata) {
    return new Promise(function(resolve,reject) {

        fs.readFile('version.json', 'utf8', function (err, data) {
            if (err) throw err; // we'll not consider error handling for now
            var json = JSON.parse(data);
            metadata['software'].push({
                'name' : 'ms-checker',
                'version' : json.revision[0],
                'run-date' : (new Date()).toISOString()
            });
            resolve();
        });
    });
};

var get_score_metadata = function(db,metadata) {
    var score_name, score_description;

    return db.each(score_type_sql,[],function(err,score) {
        score_name = score.ScoreName;
        score_description = score.Description;
    }).then(function() {
        metadata['peptide_search_engine_score'] = score_name + " : " + score_description;
    });
};

var populate_metadata = function(db) {
    var metadata = {
        'msdata-version' : MSDATA_FORMAT_VERSION 
    };

    return Promise.all( [ get_raw_filenames(db,metadata),
    get_fasta_filenames(db,metadata),
    get_self_version(metadata),
    get_pd_metadata(db,metadata),
    get_score_metadata(db,metadata) ] ).then(function() {
        return metadata;
    });
};

exports.get_metadata = populate_metadata;