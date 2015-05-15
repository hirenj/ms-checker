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
        }
    });
    grunt.loadNpmTasks('grunt-git-describe');
    grunt.loadNpmTasks('grunt-githooks');

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
}