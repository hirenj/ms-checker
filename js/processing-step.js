var util = require('util');

var ProcessingStep = function ProcessingStep() {

};

util.inherits(ProcessingStep,require('events').EventEmitter);

ProcessingStep.prototype.notify_task = function(task) {
	this.emit('task',task);
};

ProcessingStep.prototype.notify_progress = function(index,total) {
    if ( ! this.hasOwnProperty("last_frac") ) {
    	this.last_frac = 0;
    }
    var frac = parseFloat((index / (total || 1)).toFixed(2));
    if ( frac !== this.last_frac ) {
        this.emit('progress',frac);
        this.last_frac = frac;
    }
};

module.exports = exports = ProcessingStep;
