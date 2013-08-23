describe("Canvas Class for Animation", function() {
    var canvas = new Canvas(JSONObj);

    it("should be an instance of Canvas", function() {
    expect(canvas instanceof Canvas).toBeTruthy(); });


    it("should have width as in JSON object passed.", function() {
    expect(canvas.width).toEqual(JSONObj.width); })


    it("should have height as in JSON object passed.", function() {
    expect(canvas.height).toEqual(JSONObj.height); });

    it("should have a predefined Object3D for adding axes", function() {
    expect(canvas.axes instanceof THREE.Object3D).toBeTruthy(); });

    it("should have a mesh dictionary for holding shape objects", function() {
    expect(canvas.mesh_dict).toBeDefined(); });    

    // Now we initialize the canvas ...
    canvas.initialize();

    it("should have a WebGL Renderer after initialize", function() {
    expect(canvas.renderer instanceof THREE.WebGLRenderer).toBeTruthy(); });

    it("should have a scene after initialize", function() {
    expect(canvas.scene instanceof THREE.Scene).toBeTruthy(); });


    it("should have x axis after initialize", function() {
    expect(canvas.axes[0] instanceof THREE.PlaneGeometry).toBeTruthy(); });
    
    it("should have y axis after initialize", function() {
    expect(canvas.axes[0] instanceof THREE.PlaneGeometry).toBeTruthy(); });
    
    it("should have z axis after initialize", function() {
    expect(canvas.axes[0] instanceof THREE.PlaneGeometry).toBeTruthy(); });

    it("should have a Light after initialize", function() {
    expect(canvas.light instanceof THREE.PointLight).toBeTruthy(); });

    it("should have a camera defined ", function() {
    expect(canvas.camera instance of THREE.Camera).toBeTruthy(); });    

    it("should have camera controls  defined ", function() {
    expect(canvas.camera_controls instance of THREE.TrackBallControls).toBeTruthy(); });        
    
    it("should have a reset camera button defined ", function() {
    expect(canvas.reset_button).toBeDefined(); });     


    it("should have an start animation button defined." function() {
    expect(canvas.start_animation_button).toBeDefined(); });        

    it("should have a pause animation button defined." function() {
    expect(canvas.pause_animation_button).toBeDefined(); });            

    it("should have a add_frames" function() {
    expect(canvas.add_frames).toBeDefined(); });       
    
    canvas.add_frames();
    
    it("should have a add_frames" function() {
    expect(canvas.add_frames).toBeDefined(); });       
    
    s.add_frames();
    
    // Checking for mesh dictionary ..
    
    it("mesh dictionary should be populated with frames from JSON object" function() {
    expect(canvas.mesh_dict['frame1'] instanceof THREE.Mesh).toBeTruthy(); });       
    
    it("mesh dictionary member 1 should match from JSON Object." function() {
    expect(canvas.mesh_dict['frame1'] instanceof THREE.Mesh).toBeTruthy(); });       
    
    it("mesh dictionary member 2 should match from JSON Object." function() {
    expect(canvas.mesh_dict['frame2'] instanceof THREE.Mesh).toBeTruthy(); });       

    it("mesh dictionary member 3 should match from JSON Object." function() {
    expect(canvas.mesh_dict['frame3'] instanceof THREE.Mesh).toBeTruthy(); });       
    
    // Add TrackBall controls should be defined ... 
    it("should have a  method to initiate trackball controls" function() {
    expect(canvas.add_track_ball_controls).toBeDefined(); });                
        
    canvas.animation_loop = 0;
    // This should run animation once.    
    canvas.run_animation();
    

    
  });

