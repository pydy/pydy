
describe("DynamicsVisualizer's main class should  ", function() {
    it("have an initializer which returns false, when non-supported browser is used", function() {
        var browserState = Modernizr.canvas && Modernizr.webgl;
        alert(browserState);
    	var returnVal = DynamicsVisualizer.init();
        expect(returnVal).toBe(browserState);
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

describe("DynamicsVisualizer's UI ", function() {
    it("should have Play Animation button disabled, initially!", function() {
        loadFixtures("uiFixtures.html")
        expect(jQuery("#play-animation")).toBeMatchedBy('.disabled');
    });

});
