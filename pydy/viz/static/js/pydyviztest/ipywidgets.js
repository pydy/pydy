define(["jupyter-js-widgets"], function(widgets) {

    var HelloView = widgets.DOMWidgetView.extend({

        render: function() {
            this.value_changed();
            this.model.on('change:value', this.value_changed, this);
        },

        value_changed: function() {
            this.el.textContent = this.model.get('value');
        },

    });

    return {
        HelloView : HelloView
    };
});
