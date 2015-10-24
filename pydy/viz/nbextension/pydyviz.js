define([
    "widgets/js/widget",
    "underscore",
], function(widget,  _){
    console.log("[PyDy INFO] loading pydyviz");

    var TrajectoryLinkModel = widget.WidgetModel.extend({

        initialize: function() {
            this.on("change", this.update_bindings, this);
        },

        update_value: function(source, target) {
            if (this.updating) {
                return;
            }
            this.updating = true;
            try {
                if (target) {
                    var value = source.get("value");
                    var t0 = source.get("min");
                    var dt = source.get("step");
                    var index = Math.round((value - t0) / dt);
                    target.set("position", this.get("position")[index]);
                    target.set("quaternion", this.get("quaternion")[index]);
                    target.save_changes();
                }
            } finally {
                this.updating = false;
            }
        },

        update_bindings: function() {
            this.cleanup();
            this.source = this.get("source");
            this.target = this.get("target");
            if (this.source) {
                this.listenTo(this.source, "change:value", function() {
                    this.update_value(this.source, this.target);
                }, this);
                this.update_value(this.source, this.target);
                this.listenToOnce(this.source, "destroy", this.cleanup, this);
            }
            if (this.target) {
                this.listenToOnce(this.target, "destroy", this.cleanup, this);
            }
        },

        cleanup: function() {
            if (this.source) {
                this.stopListening(this.source, "change:value", null, this);
                this.stopListening(this.source, "destroy", null, this);
            }
            if (this.target) {
                this.stopListening(this.target, "destroy", null, this);
            }
        },

    }, {

        serializers: _.extend({
            target: {deserialize: widget.unpack_models},
            source: {deserialize: widget.unpack_models},
        }, widget.WidgetModel.serializers),

    });

    return {
        "TrajectoryLinkModel": TrajectoryLinkModel,
    }
});
