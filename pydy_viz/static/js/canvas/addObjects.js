Canvas.prototype.addControls = function(){
    primaryControls.update();
    renderer.render(scene, primaryCamera);
    this.controlsID = requestAnimationFrame(Canvas.prototype.addControls);
};


Canvas.prototype.resetControls = function() {
primaryControls.reset();

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
            var _Camera = new THREE.OrthographicCamera(
                                         JSONObj.width / - 2, JSONObj.width / 2, 
                                         JSONObj.height / 2, JSONObj.height / - 2,
                                         _camera.near, _camera.far );
            var _element = new Float32Array(_camera.simulation_matrix[0]);
            var initMatrix = new THREE.Matrix4();
            initMatrix.elements = _element;
            _Camera.applyMatrix(initMatrix);
        
        }
        console.log(_Camera);

        // add a small cube for a camera representation ...
        var _material = new THREE.MeshLambertMaterial(0xffffff);
        var _geom = new THREE.CubeGeometry(2,2);
        var _pointCamera = new THREE.Mesh(_geom, _material);
        _pointCamera.position = _Camera.position;
        
        this.cameraPoints.add(_pointCamera)
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
            
            // For point of light ...
            var _geom = new THREE.SphereGeometry(1,100,100);
            var _material = new THREE.MeshBasicMaterial(color);
            var _lightPoint = new THREE.Mesh(_geom, _material);
            _lightPoint.position = _Light.position
            

            
        // case "MoreLights" ....
        
        }
        
        this.lights.add(_Light);
        this.lightPoints.add(_lightPoint);

        
       
    }
    this.scene.add(this.lights);
};

Canvas.prototype.addFrames = function(){
    for(var key in JSONObj.frames){
        var _frame = JSONObj.frames[key];
        console.log(_frame);
        console.log(_frame.shape);
        var _color = new THREE.Color().setRGB(_frame.shape.color[0],
                           _frame.shape.color[1], _frame.shape.color[2]);
                                    
        var _material = new THREE.MeshLambertMaterial({
                                    color:        _frame.shape.color,
                                    opacity: 0.5
                                                     });        

        switch(_frame.shape.type){

        case "Cube":
            var _geometry = new THREE.CubeGeometry(
                                      _frame.shape.length,
                                      _frame.shape.length,
                                      _frame.shape.length,                                                                            
                                      50, 50, 50);          
            break;                          

        case "Cylinder":
            var _geometry = new THREE.CylinderGeometry(
                                      _frame.shape.radius,
                                            _frame.shape.radius,
                                      _frame.shape.length,
                                      50,50);
            break;
                                      
        case "Cone":        
            var _geometry = new THREE.CylinderGeometry(
                                      _frame.shape.radius,
                                      _frame.shape.radius/100,
                                      _frame.shape.length,
                                      50,50);        
            break;
       
        case "Sphere":        
            var _geometry = new THREE.SphereGeometry(
                                      _frame.shape.radius,
                                      100,100);
            break;

        case "Circle":        
            var _geometry = new THREE.SphereGeometry(
                                      _frame.shape.radius,
                                      8);
            break;            
           
        case "Plane":        
            var _geometry = new THREE.PlaneGeometry(
                                      _frame.shape.length,
                                      _frame.shape.width,                                      
                                      100);
            break;                        
            
        case "Tetrahedron":
            var _geometry = new THREE.TetrahedronGeometry(
                                      _frame.shape.radius);
            break;                                    
            
        case "Octahedron":
            var _geometry = new THREE.OctahedronGeometry(
                                      _frame.shape.radius);
            break;                                    
                        
        case "Icosahedron":
            var _geometry = new THREE.IcosahedronGeometry(
                                      _frame.shape.radius);
            break;                                                
            
        case "Torus":
            var _geometry = new THREE.TorusGeometry(
                                      _frame.shape.radius,
                                      _frame.shape.tube_radius,100
                                      );
            break;                                                            
            
        case "TorusKnot":
            var _geometry = new THREE.TorusKnotGeometry(
                                      _frame.shape.radius,
                                      _frame.shape.tube_radius,100
                                      );
            break;                                                                        
            
        //. ......
        }

        var _mesh = new THREE.Mesh(_geometry, _material);
        var _element = new Float32Array(_frame.simulation_matrix[0]);
        var initMatrix = new THREE.Matrix4();
        initMatrix.elements = _element;
        _mesh.matrix.identity();
        _mesh.applyMatrix(initMatrix);
        this.frames.add(_mesh);
        
       
    }
    this.scene.add(this.frames);
    for(var key in this.frames.children)
    {
    console.log("Auto Update activated");
    this.frames.children[key].matrixAutoUpdate = false;
    
    }
    
    
    
};
