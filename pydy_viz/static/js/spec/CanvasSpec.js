describe("Canvas Initializer", function() {
    var canvas = new Canvas(JSONObj);

    it("should be an instance of Canvas", function() {
    expect(canvas instanceof Canvas).toBeTruthy; }


    it("should have width  as in JSON object passed.", function() {
    expect(canvas.width).toEqual(JSONObj.width); }


    it("should have height as in JSON object passed.", function() {
    expect(canvas.height).toEqual(JSONObj.height); }

    it("should have a predefined Object3D for adding axes", function() {
    expect(canvas.axes instanceof THREE.Object3D).toBeTruthy; }

    // Now we initialize the canvas ...
    canvas.initialize();

    it("should have a WebGL Renderer after initialize", function() {
    expect(canvas.renderer instanceof THREE.WebGLRenderer).toBeTruthy; }

    it("should have a scene after initialize", function() {
    expect(canvas.scene instanceof THREE.Scene).toBeTruthy; }


    it("should have x axis after initialize", function() {
    expect(canvas.axes[0] instanceof THREE.PlaneGeometry).toBeTruthy;


  });

}
