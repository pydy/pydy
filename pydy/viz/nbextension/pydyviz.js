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

    var PlayLinkModel = widget.WidgetModel.extend({

        initialize: function() {
            this.on("change", this.update_bindings, this);
        },

        update_bindings: function() {
            this.cleanup();
            this.play = this.get("play");
            this.slider = this.get("slider");
            this.speedup = this.get("speedup");
            this.loop = this.get("loop");

            this.listenTo(this.play, "change:value", function () {
                this.on_control_click(this.play, this.slider, this.speedup, this.loop);
            }, this);
            this.play.set("value", false);
            this.play.save_changes();
        },

        on_control_click: function(control, slider, speedup, loop) {
            if (this.updating) {
                return;
            }
            this.updating = true;
            try {
                if (control && slider && speedup && loop) {
                    var go = control.get("value");
                    if (!go) {
                        window.clearTimeout(this.animationId);
                    } else {
                        var t0 = slider.get("min");
                        var tf = slider.get("max");
                        var dt = slider.get("step");
                        var end = Math.round((tf - t0) / dt);
                        this.animationId = window.setInterval(function () {
                            var time = slider.get("value");
                            // make sure we don't exceed the length of the
                            // trajectory
                            var index = Math.round((time - t0) / dt);
                            if (index >= end) {
                                slider.set("value", t0);
                                if (!loop.get("value")) {
                                    window.clearTimeout(this.animationId);
                                    control.set("value", false);
                                }
                            } else {
                                slider.set("value", time + dt);
                            }
                            slider.save_changes();
                        }, dt * speedup.get("value") * 1000);
                    }
                }
            } finally {
                this.updating = false;
            }
        },

        cleanup: function() {
            if (this.play) {
                this.stopListening(this.play, "change:value", null, this);
                this.stopListening(this.play, "destroy", null, this);
            }
            if (this.loop) {
                this.stopListening(this.loop, "destroy", null, this);
            }
            if (this.slider) {
                this.stopListening(this.slider, "destroy", null, this);
            }
            if (this.speedup) {
                this.stopListening(this.speedup, "destroy", null, this);
            }
        },

    }, {

        serializers: _.extend({
            play: {deserialize: widget.unpack_models},
            loop: {deserialize: widget.unpack_models},
            slider: {deserialize: widget.unpack_models},
            speedup: {deserialize: widget.unpack_models},
        }, widget.WidgetModel.serializers),

    });

    return {
        "TrajectoryLinkModel": TrajectoryLinkModel,
        "PlayLinkModel": PlayLinkModel,
    }
});
