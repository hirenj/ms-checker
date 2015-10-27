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
        }
    });
    grunt.loadNpmTasks('grunt-git-describe');
    grunt.loadNpmTasks('grunt-githooks');
    grunt.loadNpmTasks('grunt-node-inspector');

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