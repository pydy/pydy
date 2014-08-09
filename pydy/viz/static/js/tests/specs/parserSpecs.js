describe("DynamicsVisualizer's Parser class should have a loadScene method which when called", function() {
	beforeEach(function(done) {
        setTimeout(function() {
          DynamicsVisualizer.sceneFilePath = "sample_data/scene_desc.json";
          DynamicsVisualizer.Parser.loadScene();
          done();
        }, 1000);
        });
    it("should have DynamicsVisualizer.model defined", function(done) {
        //DynamicsVisualizer.sceneFilePath = "sample_data/scene_desc.json";
        //DynamicsVisualizer.Parser.loadScene();
        expect(DynamicsVisualizer.model).toBeDefined();
        done();
        
    });

    it("should load camera to be rendered", function() {
//        expect(DynamicsVisualizer._scene[0]).toBe(THREE.PerspectiveCamera);
    });

    it("should load light to be rendered", function() {
        
        //expect(DynamicsVisualizer._scene[1]).toBe(THREE.PointLight});
    });

    it("should load axes to be rendered", function() {
        //expect(DynamicsVisualizer._scene[2]).toBe(THREE.AxisHelper);
    });

    it("play animation should be un-disabled", function() {
//        expect(jQuery("#play-animation")).not.toBeMatchedBy('.disabled');
    });
    it("Show Model should be un-disabled", function() {
//        expect(jQuery("#show-model")).not.toBeMatchedBy('.disabled');
    });

    it("should save simulation data in an object upon success", function() {
//        expect(DynamicsVisualizer.simData).toBeDefined();
        
    });

    it("should create a timeArray on succesful completion", function() {
//        expect(DynamicsVisualizer._timeArray).toBeDefined();
    });

});
