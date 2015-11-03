var fs = require('fs');


var contaminants = { 'empty' : true, 'identifiers' : [], 'version' : 'NA', 'retrieved' : 'NA' };

var read_contaminants = function(done) {
    if ( ! contaminants.empty ) {
        return done();
    }
    fs.readFile('contaminants.json', 'utf8', function (err, data) {
        if (err) throw err; // we'll not consider error handling for now
        var json = JSON.parse(data);
        contaminants = json;
        exports.identifiers = contaminants.identifiers;
        done();
    });
}

var get_version = function(metadata) {
    return ( new Promise(function(resolve,reject) {
        read_contaminants(function() {
            resolve();
        });
    })).then(function() {
        metadata['contaminants'] = { 'source' : 'maxquant.org', 'updated' : contaminants.version, 'retrieved' : contaminants.retrieved };
    });
};

read_contaminants(function() {});

exports.get_version = get_version;
exports.identifiers = contaminants.identifiers;