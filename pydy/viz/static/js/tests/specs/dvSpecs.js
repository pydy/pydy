
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
TODO add ui fixtures, and test UI
describe("DynamicsVisualizer's UI ", function() {
    beforeEach(function(){
       loadFixtures("sample_data/uiFixtures.html")
    });
    it("should have Play Animation button disabled, initially!", function() {
        expect(jQuery("#play-animation")).toBeMatchedBy('.disabled');
    });

});

**/

