var peptide_search = require('./peptide.js');
var util = require('util');

var Ambiguous = function Ambiguous() {

};

util.inherits(Ambiguous,require('./processing-step.js'));

module.exports = exports = new Ambiguous();


const retrieve_ambiguous_peptides_glyco_sql = 'SELECT \
    Peptides.PeptideID, \
    Peptides.SpectrumID, \
    Peptides.Sequence, \
    Peptides.SearchEngineRank, \
    SpectrumHeaders.ScanNumbers as scan, \
    SpectrumHeaders.Mass as mass, \
    SpectrumHeaders.RetentionTime as retentionTime, \
    SpectrumHeaders.Charge as charge, \
    FileInfos.Filename, \
    ScanEvents.ActivationType as ActivationType \
FROM FileInfos \
    LEFT JOIN MassPeaks USING(FileID) \
    LEFT JOIN SpectrumHeaders USING(MassPeakID) \
    LEFT JOIN ScanEvents USING(ScanEventID) \
    JOIN Peptides using(SpectrumID) \
    WHERE \
        (   FileInfos.FileName LIKE "%xGalNAc%" \
        OR \
            FileInfos.FileName LIKE "%xHexNAc%" \
        OR \
            FileInfos.FileName LIKE "%Hex%" \
        OR \
            FileInfos.FileName LIKE "%Man%" \
        OR \
            FileInfos.FileName LIKE "%xT.%" \
        OR \
            FileInfos.FileName LIKE "%xTn.%" \
        ) AND \
            FileInfos.Filename NOT LIKE "%raw" \
    AND \
        Peptides.ConfidenceLevel = 3 \
    AND \
        Peptides.SearchEngineRank <= 4 \
';

const retrieve_ambiguous_peptides_sql = 'SELECT \
    Peptides.PeptideID, \
    Peptides.SpectrumID, \
    Peptides.Sequence, \
    Peptides.SearchEngineRank, \
    SpectrumHeaders.ScanNumbers as scan, \
    SpectrumHeaders.Mass as mass, \
    SpectrumHeaders.RetentionTime as retentionTime, \
    SpectrumHeaders.Charge as charge, \
    FileInfos.Filename \
FROM FileInfos \
    LEFT JOIN MassPeaks USING(FileID) \
    LEFT JOIN SpectrumHeaders USING(MassPeakID) \
    JOIN Peptides using(SpectrumID) \
    WHERE \
        Peptides.ConfidenceLevel = 3 \
    AND \
        Peptides.SearchEngineRank <= 4 \
';


var clean_comp = function(comp) {
    comp = comp.replace('xTn','xHexNAc');
    comp = comp.replace('xT','xHexHexNAc');
    comp = comp.replace('xGalNAc','xHexNAc');
    return comp;
}

var retrieve_ambiguous_peptides = function(db) {
    var self = this;
    var sql = this.conf.get("match-all-ambiguous") ? retrieve_ambiguous_peptides_sql : retrieve_ambiguous_peptides_glyco_sql;

    self.notify_task('Finding ambiguous peptides');
    self.notify_progress(0,1);
    return db.all(sql).then(function(peps) {
        self.notify_progress(1,1);
        self.notify_task('Populating peptide metadata');
        self.notify_progress(0,peps.length);

        var idx = 0;
        var update_fractions = function(pep) {
            idx += 1;
            self.notify_progress(idx,peps.length);
        };

        total_peps = peps.length;
        return Promise.all( peps.length > 0 ? peps.map(function(pep) { return peptide_search.produce_peptide_data(db,pep).then(update_fractions); }) : [true] ).then(function() {
            self.notify_progress('progress',peps.length,peps.length);
            peptide_search.cleanup(db);
            peps.forEach(function(pep) {
                var composition = (pep.FileName || "").match(/\d+\x[^%]+(?=\.mgf)/);
                composition = Array.prototype.concat.apply( [], composition.map(function(comp) { return comp.split(/[_,]/); }) );
                if ( composition ) {
                    pep.Composition = composition.map(function(comp) { return comp.toString(); }).map(clean_comp);
                }
                delete pep.FileName;
            });
            console.log("Total "+peps.length+" ambiguous peptides");
            return peps;
        });
    } );
};

exports.retrieve_peptides = retrieve_ambiguous_peptides;