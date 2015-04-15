/* These are the methods, which cannot be tested on console
 * Due to absence of a renderer in the bash:
 * 
 *  - DynamicsVisualizer.Scene._createRendere();
 *  - DynamicsVisualizer.Scene._addTrackBallControls();
*/ 

describe("DynamicsVisualizer's Scene class should  ", function() {

	beforeEach(function(){
        DynamicsVisualizer.sceneFilePath = "sample_data/scene_desc.json";
        DynamicsVisualizer.Scene._createEmptyScene();
    });
    it("have a method which could initialize Scene", function() {
        expect(DynamicsVisualizer._scene).toBeDefined();
    });
    it("have a method which could initialize Axes", function() {
    	DynamicsVisualizer._addAxes();
        expect(DynamicsVisualizer._scene.children[0]).toEqual(jasmine.any(THREE.AxisHelper));;
    });
    it("have a method which could initialize Default camera", function() {
    	DynamicsVisualizer._addDefaultCamera();
        expect(DynamicsVisualizer._scene.children[0]).toEqual(jasmine.any(THREE.PerspectiveCamera));
    });

    it("have a method which could initialize Default Light", function() {
    	DynamicsVisualizer._addDefaultLight();
        expect(DynamicsVisualizer._scene.children[0]).toEqual(jasmine.any(THREE.PointLight));;
    });
    


});
