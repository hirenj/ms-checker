(function(scope) {
    var promisify_sqlite = function(db) {
        var old_all = db.all;
        db.all = function(sql,vals) {
            var args = Array.prototype.splice.call(arguments);
            return new Promise(function(resolve,reject) {
                old_all.call(db,sql,vals,function(err,out_vals) {
                    if (err) {
                        reject(err);
                    } else {
                        resolve(out_vals);
                    }
                });
            });
        };
        var old_each = db.each;
        db.each = function(sql,vals,callback) {
            if (! callback ) {
                callback = function() {};
            }
            var args = Array.prototype.splice.call(arguments);
            return new Promise(function(resolve,reject) {
                old_each.call(db,sql,vals,callback,function(err,vals) {
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
            if ( ! cached_statements[sql] ) {
                cached_statements[sql] = db.prepare(sql);
            }
            cached_statements[sql].active_count = (cached_statements[sql].active_count || 0) + 1;
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
            if (cached_statements[sql] && cached_statements[sql].active_count > 1) {
                cached_statements[sql].active_count -= 1;
                return;
            }
            if (cached_statements[sql]) {
                cached_statements[sql].finalize();
            }
            delete cached_statements[sql];
        };
    };
    module.exports = promisify_sqlite;
})(module.exports);