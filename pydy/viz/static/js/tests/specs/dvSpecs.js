describe("DynamicsVisualizer's main class should  ", function() {
    it("have an initializer which returns false, when non-supported browser is used", function() {
        var browserState = Modernizr.canvas && Modernizr.webgl;
        alert(browserState);
    	var returnVal = DynamicsVisualizer._initialize();
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
    beforeEach(function() {
        loadFixtures("uiFixtures.html");
    });

    it("should have Play Animation button disabled, initially!", function() {
        expect(jQuery("#play-animation")).toBeMatchedBy('.disabled');
    });

    it("should have Play Animation button disabled, initially!", function() {
        expect(jQuery("#play-animation")).toBeMatchedBy('.disabled');
    });

    it("should have Stop Animation button hidden, initially!", function() {
        expect(jQuery("#stop-animation")).toBeHidden();
    });
});

describe("After activating UI controls,", function() {
    beforeEach(function() {
        loadFixtures("uiFixtures.html");
        DynamicsVisualizer._initialize();
        DynamicsVisualizer.activateUIControls();
    });

    it("On clicking Load Simulation button, its click function should be triggered", function() {
        var spyEvent = spyOnEvent('#simulation-load', 'click')
        jQuery("#simulation-load").click()
        expect('click').toHaveBeenTriggeredOn('#simulation-load');
    });
    it("On clicking Load Simulation button, Scene File path should be updated", function() {
        jQuery("#simulation-load").click();
        expect(DynamicsVisualizer.sceneFilePath).toEqual('sample_data/scene_desc.json');
    });
});
