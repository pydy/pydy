
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
    
});

describe("DynamicsVisualizer's main class should  ", function() {
    beforeEach(function(){
        DynamicsVisualizer.sceneFilePath = "sample_data/scene_desc.json";
    });
    
    it("have a method to return the basepath of the simulation file url", function() {
        var expectedBasePath = "sample_data/";
        expect(DynamicsVisualizer.getBasePath()).toBe(expectedBasePath);
    });

    it("have a method to return the file extension of the simulation file url", function() {
        var expectedFileExtension = "json";
        expect(DynamicsVisualizer.getFileExtenstion()).toBe(expectedFileExtension);
    });
    
});

/**
Bug in jasmine-jquery, added an issue at:
https://github.com/velesin/jasmine-jquery/issues/199

//TODO add ui fixtures, and test UI
describe("DynamicsVisualizer's UI ", function() {
    beforeEach(function () {
    jQuery.ajax({
      async: false, // must be synchronous to guarantee that no tests are run before fixture is loaded
      dataType: 'html',
      url: 'sample_data/uiFixture.html',
      success: function(data) {
        $('body').append($(data));
      }
    });
  });
    
    it("should have Play Animation button disabled, initially!", function() {
        loadFixtures("uiFixtures.html")
        expect(jQuery("#play-animation")).toBeMatchedBy('.disabled');
    });

});

**/