var fs = require('fs');
var nconf = require('nconf');

var contaminants = require('./contaminants');

const fasta_files_sql = "SELECT FileName FROM FastaFiles";

const raw_files_sql = "SELECT FileName, FileTime FROM FileInfos";

const pd_metadata_sql = "SELECT SoftwareVersion, Date FROM SchemaInfo";

const score_type_sql = 'SELECT ScoreID, ScoreName, Description FROM ProcessingNodeScores';

const MSDATA_FORMAT_VERSION = "1.4";

var get_raw_filenames = function(db,metadata){
    if ( ! db ) {
        return Promise.resolve();
    }
    var raw_filenames = [];
    return db.each(raw_files_sql,[],function(err,file) {
        raw_filenames.push({'file': file.FileName, 'created': new Date(file.FileTime).toISOString() });
    }).then(function() {
        metadata['ms_run-location'] = raw_filenames;
    });
};

var get_fasta_filenames = function(db,metadata) {
    if ( ! db ) {
        return Promise.resolve();
    }
    var fasta_filenames = [];
    return db.each(fasta_files_sql,[],function(err,file) {
        fasta_filenames.push(file.FileName);
    }).then(function() {
        metadata['database'] = fasta_filenames;
        metadata['database_version'] = fasta_filenames.map(function() {  return "" });
    });
};

var get_sample_metadata = function(metadata) {
    metadata['sample'] = {};
    metadata['sample']['species'] = nconf.get('source-organism');
    metadata['sample']['tissue'] = nconf.get('source-organism_part');
    if (nconf.get('source-cell_line')) {
        metadata['sample']['cell_type'] = nconf.get('source-cell_line');
        metadata['sample']['cell_type_id'] = 'RRID:'+nconf.get('source-cellosaurus_id');
    }
    if (nconf.get('perturbation-wt')) {
        metadata['sample']['wt'] = nconf.get('perturbation-wt').map((ko) => "entrez:"+ko.entrez );
    }
    if ((nconf.get('perturbation-ko') || '').length) {
        metadata['sample']['ko'] = nconf.get('perturbation-ko').map((ko) => "entrez:"+ko.entrez );
    }
    if ((nconf.get('perturbation-ki') || '').length) {
        metadata['sample']['ki'] = nconf.get('perturbation-ki').map((ki) => "entrez:"+ki.entrez );
    }
    if ((nconf.get('perturbation-other') || '').length) {
        metadata['sample']['perturbation-other'] = nconf.get('perturbation-other');
    }
    metadata['channel_samples'] = {};
    for (let channel_sample of nconf.get('perturbations') || []) {
        metadata['channel_samples'][channel_sample['perturbation-identifier']] = Object.assign({},channel_sample);
        delete metadata['channel_samples'][channel_sample['perturbation-identifier']]['perturbation-identifier'];
    }
};

var get_pd_metadata = function(db,metadata) {
    if ( ! db ) {
        return Promise.resolve();
    }
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

var get_ppm_cutoffs = function(metadata) {
    var conf = nconf.get();
    metadata['ppm-cuttofs'] = {};
    Object.keys(conf).filter((key) => key.indexOf('ppm') == 0).forEach(function(key) {
        metadata['ppm-cuttofs'][key] = conf[key];
    });
};


var write_feature_flags = function() {
    var opts = nconf.get();
    var results = [];
    Object.keys(opts).forEach(function(flag) {
        if (flag.match('feature')) {
            results.push(flag);
        }
    });
    return results.join(',');
};

var get_quant_conf = function(metadata) {
    var opts = nconf.get();
    var results = { };
    Object.keys(opts).forEach(function(flag) {
        if (flag.match('quant-channel-')) {
            results[flag.split('-')[2]] = opts[flag];
        }
    });
    if (Object.keys(results).length > 0) {
        metadata.quantitation = { 'channels' : results };
    }
};


var get_self_version = function(metadata) {
    return new Promise(function(resolve,reject) {

        fs.readFile(__dirname+'/../version.json', 'utf8', function (err, data) {
            if (err) throw err; // we'll not consider error handling for now
            var json = JSON.parse(data);
            metadata['software'] = metadata['software'] || [];
            metadata['software'].push({
                'name' : 'ms-checker',
                'version' : json.revision[0],
                'run-date' : (new Date()).toISOString(),
                'feature-flags' : write_feature_flags()
            });
            resolve();
        });
    });
};

var get_doi = function(metadata) {
    if (nconf.get('doi')) {
        metadata.doi = nconf.get('doi');
    }
};

var get_score_metadata = function(db,metadata) {
    if ( ! db ) {
        return Promise.resolve();
    }
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
        'msdata-version' : MSDATA_FORMAT_VERSION,
        'mimetype' : 'application/json+msdata',
        'tags' : 'ccg'
    };

    return Promise.all( [ get_raw_filenames(db,metadata),
    get_fasta_filenames(db,metadata),
    get_self_version(metadata),
    get_pd_metadata(db,metadata),
    get_sample_metadata(metadata),
    get_ppm_cutoffs(metadata),
    get_doi(metadata),
    get_quant_conf(metadata),
    contaminants.get_version(metadata),
    get_score_metadata(db,metadata) ] ).then(function() {
        return metadata;
    });
};

var combine_metadata = function(metadatas) {
    return metadatas.reduce(function(result,next) {
        if ( ! result ) {
            return next;
        }
        if (next['ms_run-location']) {
            result['ms_run-location'] = (result['ms_run-location'] || []).concat(next['ms_run-location']);
        }
        return result;
    });
};


exports.version = get_self_version;

exports.get_metadata = populate_metadata;

exports.combine = combine_metadata;
