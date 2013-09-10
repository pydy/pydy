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
            var _color = new THREE.Color().setRGB(_frame.shape.color[0],
                           _frame.shape.color[1], _frame.shape.color[2])
                                    
            var _material = new THREE.MeshLambertMaterial({
                                    color:        _frame.shape.color,
                                    wireframe:          true,
                                    wireframeLinewidth: 0.1,
                                    opacity: 0.5
                                                     })
                                                                    
            var _geometry = new THREE.CylinderGeometry(
                                      _frame.shape.radius,
                                      _frame.shape.radius,
                                      _frame.shape.height,
                                      50,50);
            var _mesh = new THREE.Mesh(_geometry, _material);
            
            var _element = new Float32Array(_frame.simulation_matrix[0]);
            var initMatrix = new THREE.Matrix4();
            initMatrix.elements = _element;
            _mesh.matrix.identity();
            _mesh.applyMatrix(initMatrix);
            // case for more shapes ...
        
        }
        this.frames.add(_mesh);
        
       
    }
    this.scene.add(this.frames);
    for(var key in this.frames.children)
    {
    console.log("Auto Update activated");
    this.frames.children[key].matrixAutoUpdate = false;
    
    }
    
    
    
};
