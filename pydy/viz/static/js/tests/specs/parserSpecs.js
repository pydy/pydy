describe("DynamicsVisualizer's Parser class should have a loadScene method which when called", function() {
	beforeEach(function(){
        DynamicsVisualizer.sceneFilePath = "scene_desc.json";
        DynamicsVisualizer.Parser.loadScene();
    });
    
    
    it("should load objects to be rendered", function() {
        console.log(DynamicsVisualizer.model);
        //expect(DynamicsVisualizer._scene).toBe(expectedBasePath);
    });

    it("should load lights to be rendered", function() {
        
        //expect(DynamicsVisualizer._scene).toBe(expectedBasePath);
    });

    it("should load cameras to be rendered", function() {
        
        //expect(DynamicsVisualizer._scene).toBe(expectedBasePath);
    });

});

describe("On succesful calling Parser.loadscene, loadUIElements should be called and", function() {
	beforeEach(function(){
        DynamicsVisualizer.sceneFilePath = "sample_data/scene_desc.json";
        DynamicsVisualizer.Parser.loadScene();
    });
    
    
    it("play animation should be un-disabled", function() {
        expect(jQuery("#play-animation")).not.toBeMatchedBy('.disabled');
    });
    it("Show Model should be un-disabled", function() {
        expect(jQuery("#show-model")).not.toBeMatchedBy('.disabled');
    });

    it("Code mirror's editor object should be defined", function() {
        expect(DynamicsVisualizer.editor).toBeDefined();
    });
});




describe("DynamicsVisualizer's Parser class should have a loadSimulation method which ", function() {
	beforeEach(function(){
        DynamicsVisualizer.sceneFilePath = "sample_data/scene_desc.json";
        DynamicsVisualizer.loadScene();
        DynamicsVisualizer.Parser.loadSimulation();
    });
    
    it("should save simulation data in an object upon success", function() {
        
        //expect(DynamicsVisualizer._scene).toBe(expectedBasePath);
    });

});

describe("DynamicsVisualizer's Parser class should have a loadSimulation method which ", function() {
	beforeEach(function(){
        DynamicsVisualizer.sceneFilePath = "sample_data/scene_desc.json";
        DynamicsVisualizer.Parser.loadScene();
    });
    
    it("should create a timeArray on succesful completion", function() {
        
        //expect(DynamicsVisualizer._scene).toBe(expectedBasePath);
    });

});