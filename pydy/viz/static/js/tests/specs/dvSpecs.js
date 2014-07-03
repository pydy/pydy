
describe("DynamicsVisualizer's main class should  ", function() {
    
  
    it("have a webgl checker method, to check if browser is compatible", function() {
    	var browserState = Modernizr.canvas && Modernizr.webgl;
    	var returnVal = DynamicsVisualizer.isWebGLCompatible();
        expect(returnVal).toBe(browserState);
    });

    it("have an initializer which returns false, when non-supported browser is used", function() {
        var browserState = Modernizr.canvas && Modernizr.webgl;
    	var returnVal = DynamicsVisualizer.init();
        expect(returnVal).toBe(false);
    });
    /*
    TODO: Figure this out here!
    it("have a method to activate click functions for UI buttons", function() {
    	var spyEvent = spyOnEvent('#play-animation', 'show')
        jQuery('#play-animation').show();
        expect('show').toHaveBeenTriggeredOn('#play-animation');
        expect(spyEvent).toHaveBeenTriggered()
    });
    */





});


describe("DynamicsVisualizer's Parser class should  ", function() {
    
  
});


describe("DynamicsVisualizer's Scene class should  ", function() {
  

});

describe("DynamicsVisualizer's ParamEditor class should  ", function() {
    
  

});



