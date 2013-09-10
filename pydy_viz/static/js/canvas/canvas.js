// Wont work until camera, lights etc. position object is defined
// in python. ..

var Canvas = function(JSONObj) {

    this.grid = new THREE.Object3D();
    this.grid.name = "Grid";
    
    this.cameras = new THREE.Object3D()
    this.cameras.name = "Cameras";
    
    this.lights = new THREE.Object3D()
    this.lights.name = "Lights";
    
    this.frames = new THREE.Object3D();
    this.frames.name = "Frames";
    
    this.animationCounter = 0;
    var renderer, scene, primaryCamera, primaryControls

};
    
Canvas.prototype.initialize = function(){



	this.renderer = new THREE.WebGLRenderer();
	this.renderer.setSize(JSONObj.width, JSONObj.height);
	
	this.container = $('#container');
	this.container.append(this.renderer.domElement);	
	
	this.scene = new THREE.Scene();
	
	var axesMaterial = new THREE.MeshLambertMaterial(
	                            {
	                                color: 0xFFFFFF
	                            });

    var _pcamera = JSONObj.cameras[0]
    this.aspect_ratio = JSONObj.width / JSONObj.height;
    
    switch(_pcamera.type){
    
        case "PerspectiveCamera":
            this.primaryCamera = new THREE.PerspectiveCamera(
                                           _pcamera.fov, this.aspect_ratio, 
                                           _pcamera.near, _pcamera.far );
            this.primaryCamera.matrix.elements = _pcamera.simulation_matrix[0];
            this.primaryCamera.position.x = 00 //_pcamera.position[0];
            this.primaryCamera.position.y = 00 //_pcamera.position[1];
            this.primaryCamera.position.z = 100 //_pcamera.position[2];
            
    }

    // setting up primary trackball controls ...
    this.primaryControls = new THREE.TrackballControls(
                                            this.primaryCamera,
                                            this.renderer.domElement);
    
    var gridYZ = new THREE.GridHelper(100, 5);
	gridYZ.position.set( 0,0,0 );
    gridYZ.material.color = new THREE.Color(0xFFFFFF);
	gridYZ.rotation.y = Math.PI/2;
	
    this.grid.add(gridYZ);
    this.scene.add(this.grid);
    
    this.axes = new THREE.AxisHelper(JSONObj.height);
    this.axes.name = "Axes";
    this.scene.add(this.axes);
    this.renderer.render(this.scene, this.primaryCamera);
    
    // copying to vars for 
    
    primaryControls = this.primaryControls; 
    primaryCamera = this.primaryCamera;
    renderer = this.renderer;
    scene = this.scene;
};

Canvas.prototype.addControls = function(){
    primaryControls.update();
    renderer.render(scene, primaryCamera);
    this.animationID = requestAnimationFrame(Canvas.prototype.addControls);
};



Canvas.prototype.addCameras = function() {
    for(var key in JSONObj.cameras){
        var _camera = JSONObj.cameras[key];
        
        switch(_camera.type){
        
        case "PerspectiveCamera":
            var _Camera = new THREE.PerspectiveCamera(
                                       _camera.fov, this.aspect_ratio,
                                         _camera.near, _camera.far );
            var _element = new Float32Array(_camera.simulation_matrix[0]);
            var initMatrix = new THREE.Matrix4();
            initMatrix.elements = _element;
            _Camera.applyMatrix(initMatrix);
        case "OrthoGraphicCamera":
            //TODO
        
        }
        console.log(_Camera);
        this.cameras.add(_Camera);
       
    }
    this.scene.add(this.cameras);
};

Canvas.prototype.addLights = function() {
    for(var key in JSONObj.lights){
        var _light = JSONObj.lights[key];
        
        switch(_light.type){
        
        case "PointLight":
            var color = new THREE.Color().setRGB(_light.color[0],
                                    _light.color[1], _light.color[2])
            var _Light = new THREE.PointLight(color);
            var _element = new Float32Array(_light.simulation_matrix[0]);
            var initMatrix = new THREE.Matrix4();
            initMatrix.elements = _element;
            _Light.applyMatrix(initMatrix);
            
        // case "MoreLights" ....
        
        }
        this.lights.add(_Light);
       
    }
    this.scene.add(this.lights);
};

Canvas.prototype.addFrames = function(){
    for(var key in JSONObj.frames){
        var _frame = JSONObj.frames[key];
        console.log(_frame);
        console.log(_frame.shape);
        switch(_frame.shape.type){
        
        case "Cylinder":
            var color = new THREE.Color().setRGB(_frame.shape.color[0],
                           _frame.shape.color[1], _frame.shape.color[2])
                                    
            var _material = new THREE.MeshLambertMaterial({
                                     color:       _frame.shape.color,
                                     wireframe:   true,
                                     wireframeLinewidth: 0.1,
                                     opacity: 0.5
                                     });

            var _geometry = new THREE.CylinderGeometry(
                                      _frame.shape.radius,
                                      _frame.shape.radius,
                                      _frame.shape.height,
                                      50,50);
            var _mesh = new THREE.Mesh(_geometry, _material);
            
            var _element = new Float32Array(_frame.simulation_matrix[0]);
            var initMatrix = new THREE.Matrix4();
            initMatrix.elements = _element;
            _mesh.applyMatrix(initMatrix);
            
        // case for more shapes ...
        
        }
        this.frames.add(_mesh);
       
    }
    this.scene.add(this.frames);
};



Canvas.prototype.startAnimation = function(){
    for(var i = 0; i< 100; i++){
    for(var key in this.frames.children){
        
        var _element = new Float32Array(
               JSONObj.frames[key].simulation_matrix[i]
                    );
        var _matrix = new THREE.Matrix4();
        _matrix.elements = _element;
        this.frames.children[key].applyMatrix(_matrix)
        i++;
        if(i >=100) {i = 0; }
        
        }
    }
};
