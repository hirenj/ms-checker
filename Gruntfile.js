module.exports = function(grunt) {
    grunt.initConfig({
        pkg: grunt.file.readJSON('package.json'),

        'git-describe': {
            options: {
            },
            me : {

            }
        },
        githooks: {
            all: {
                'post-commit': 'version'
            }
        },
        'node-inspector': {
            dev: {}
        },
        curl: {
            contaminants: {
                src: 'http://lotus1.gwdg.de/mpg/mmbc/maxquant_input.nsf/7994124a4298328fc125748d0048fee2/$FILE/contaminants.fasta',
                dest: 'tmp/contaminants.fasta'
            }
        },
        unzip: {
            contaminants: {
                src: 'tmp/contaminants.zip',
                dest: 'tmp/'
            }
        }
    });
    grunt.loadNpmTasks('grunt-git-describe');
    grunt.loadNpmTasks('grunt-githooks');
    grunt.loadNpmTasks('grunt-node-inspector');
    grunt.loadNpmTasks('grunt-curl');
    grunt.loadNpmTasks('grunt-zip');

    grunt.registerTask('saveRevision', function() {
        grunt.event.once('git-describe', function (rev) {
            grunt.log.writeln("Git Revision: " + rev);
            grunt.option('gitRevision', rev);
        });
        grunt.task.run('git-describe');
    });

    grunt.registerTask('tag-revision', 'Tag the current build revision', function () {
        grunt.task.requires('git-describe');

        grunt.file.write('version.json', JSON.stringify({
            version: grunt.config('pkg.version'),
            revision: grunt.option('gitRevision'),
            date: grunt.template.today()
        }));
    });

    grunt.registerTask('download-contaminants', 'Download the contaminants list',function() {
        grunt.event.once('curl',function(response) {
            grunt.config.set('contaminants.updated', new Date(response.headers['last-modified']).toISOString().slice(0,10));
            grunt.config.set('contaminants.retrieved', new Date(response.headers['date']).toISOString().slice(0,10));
        });
        grunt.task.run('curl:contaminants');
    });

    grunt.registerTask('parse-contaminants',"Parse contaminants list", function() {
        grunt.task.requires('download-contaminants');
        var file = grunt.file.read('tmp/contaminants.fasta');
        var identifiers = file.split("\n").filter(function(line) { return line.match(/^>/); }).map(function(line) {
            return line.split(/\s/)[0].substring(1);
        });
        var result = { "version" : grunt.config.get('contaminants.updated'), "retrieved" : grunt.config.get('contaminants.retrieved'), "identifiers" : identifiers };
        grunt.file.write('contaminants.json',JSON.stringify(result));
    });

    grunt.registerTask('update-contaminants','Update contaminants',['download-contaminants','parse-contaminants']);

    grunt.registerTask('version', ['saveRevision', 'tag-revision']);
    grunt.registerTask('debug', function() {
        var done = this.async();
        grunt.util.spawn({
            cmd: 'node',
            args: ['--debug-brk', 'js/main.js', '--output blah', '--source blah', '"/Scratch/MSF-test/CHO T1KOT2KO/150714_Velos_(LC1)_Yoshiki_CHO_T1_T2KO_sec_VVA(june2016)_2h_HCD_ETD_IEF1_12(ETD).msf"'],
            opts: {
                //cwd: current workin directory
            }
        },
        function (error, result, code) {
            if (error) {
                grunt.log.write (result);
                grunt.fail.fatal(error);
            }
            done();
        });
        grunt.log.writeln ('node started');
    });
}