
DynamicsVisualizer._initialize();
DynamicsVisualizer.activateUIControls();
DynamicsVisualizer.Scene.create();

//  Activate TrackBall Controls..
function activateTrackballControls(){
    DynamicsVisualizer.primaryControls.update()
    DynamicsVisualizer.Scene.webgl_renderer.render(DynamicsVisualizer.Scene._scene, DynamicsVisualizer.currentCamera);
    requestAnimationFrame(activateTrackballControls);
};
activateTrackballControls();

// for bootstrap proto panga..
(function() {
    var isBootstrapEvent = false;
    if (window.jQuery) {
        var all = jQuery('*');
        jQuery.each(['hide.bs.dropdown',
            'hide.bs.collapse',
            'hide.bs.modal',
            'hide.bs.tooltip'], function(index, eventName) {
            all.on(eventName, function( event ) {
                isBootstrapEvent = true;
            });
        });
    }
    var originalHide = Element.hide;
    Element.addMethods({
        hide: function(element) {
            if(isBootstrapEvent) {
                isBootstrapEvent = false;
                return element;
            }
            return originalHide(element);
        }
    });
})();
