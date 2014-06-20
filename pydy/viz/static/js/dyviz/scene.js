

DynamicsVisualizer.Scene = Object.extend(DynamicsVisualizer, {
	
	create: function(){

		/** 
		  * This method creates the scene from the self.model
		  * and renders it onto the canvas.
		  * 
		**/ 
		var self = this;
        self._createRenderer();
        self._addDefaultLightsandCameras();
        self._addAxes();
        self._addTrackBallControls(); // and render too.
        
    
	},

	_createRenderer: function(){
		/**
		  * Creates a webGL Renderer
		  * with a default background color.
		  *
		**/ 
		var self = this;
		self.webgl_renderer = new THREE.WebGLRenderer();
        self.webgl_renderer.setSize(640, 480);
        
	    var backgroundColor = new THREE.Color(161192855); // WhiteSmoke
	    self.webgl_renderer.setClearColor(backgroundColor);	
	    var container = jQuery('#renderer');
	    container.append(self.webgl_renderer.domElement);	
        
        // new Scene..
        self._scene  = new THREE.Scene();

    },    
        
    
	_addDefaultLightsandCameras: function(){
        // This is an auxillary function 
        // for development purposes..
        // should be removed after everything 
        // is in place and working.
		var self = this;
        self._cameras = {};

		var camera = new THREE.PerspectiveCamera();
        camera.position.x = 00;
        camera.position.y = 00;
        camera.position.z = 100;
        self._cameras["init_camera"] = camera;
        self._scene.add(self._cameras.init_camera);

        self._lights = {};
        var light = new THREE.PointLight(0xffffff);
        light.position.set(10,10,-10);
        self._lights["init_light"] = light;
        self._scene.add(self._lights.init_light);

	},
    

	_addAxes: function(){
        var self = this;
        self._meshes = {};
		var self = this;
		var axes = new THREE.AxisHelper(100);
        self._meshes["axes"] = axes;
        self._scene.add(self._meshes["axes"]);
        
	},

    _addTrackBallControls: function(){
    	
    	this.primaryControls = new THREE.TrackballControls(this._cameras.init_camera,
                                            this.webgl_renderer.domElement);
    
    },

    _resetControls: function(){
    	this.primaryControls.reset();
    },

    addObjects: function(){
        var self = this;
        var objects = self.model.objects;
        for(var i in objects)  self._addIndividualObject(objects[i]);
        
    },

    addCameras: function(){
        var self = this;
        var cameras = this.model.cameras;
        for(var i in cameras)  self._addIndividualCamera(cameras[i]);

        //for(var i in self._cameras)  self._scene.add(self._cameras[i]);
    },

    addLights: function(){
        var self = this;
        var lights = this.model.lights;
        for(var i in lights)  self._addIndividualLight(lights[i]);

        //for(var i in self._lights)  self._scene.add(self._lights[i]);
    },

    _addIndividualObject: function(object){
        var self = this;
        var type = object.type;
        console.log("In add-Ind");
        console.log(object);
        
        // This is for shapes,
        // Need to do something else
        // for camera, lights
        //var material = self.Materials[object.material];
        // meanwhile..
        var material = new THREE.MeshLambertMaterial();
        switch(type) {

            case "Mesh":
                //TODO
                break;

            case "Cube":
                var geometry = new THREE.CubeGeometry(
                                  object.length,
                                  object.length,
                                  object.length, 
                                  50, 50, 50);
                break;

            case "Sphere":
                var geometry = new THREE.SphereGeometry(
                                               object.radius, 8);
                break;

            case "Cylinder":
                var geometry = new THREE.CylinderGeometry(object.radius,
                                                          object.radius,
                                                          object.length);

                break;

            // TODO: Add rest of objects too.    
        }

        var mesh = new THREE.Mesh(geometry, material);
        var element = new Float32Array(object.init_orientation);
        var initMatrix = new THREE.Matrix4();
        initMatrix.elements = element;
        mesh.matrix.identity();
        mesh.applyMatrix(initMatrix);
        mesh["object-info"] = object;
        mesh.name = object.simulation_id;
        self._scene.add(mesh);
        //self._meshes[object.simulation_id] = mesh;

        // This info is for object editing dialog..
        //self._meshes[object.simulation_id]["object-info"] = object;
    },

    _addIndividualCamera: function(camera){
        
        var self = this;
        switch(camera.type){
            case "PerspectiveCamera":
                var _camera = new THREE.PerspectiveCamera(camera.fov, 1,
                                                 camera.near, camera.far);
                var element = new Float32Array(camera.init_orientation);
                var initMatrix = new THREE.Matrix4();
                initMatrix.elements = element;
                _camera.applyMatrix(initMatrix);
                break;
            
            /*case "OrthoGraphicCamera":
            // TODO adjust JSONObj.width and height here..
                var _camera = new THREE.OrthographicCamera(
                                     JSONObj.width / - 2, JSONObj.width / 2, 
                                     JSONObj.height / 2, JSONObj.height / - 2,
                                     _camera.near, _camera.far );
                var _element = new Float32Array(_camera.simulation_matrix[0]);
                var initMatrix = new THREE.Matrix4();
                initMatrix.elements = _element;
                _camera.applyMatrix(initMatrix);
            */    
        }
        _camera.name = camera.simulation_id;
        _camera["object-info"] = camera;
        self._scene.add(_camera);

    },

    _addIndividualLight: function(light){
        // TODO: add individual lights to self._lights
        var self = this;
        var type = light.type;
        switch(light.type) { 

            case "PointLight":
                var color = new THREE.Color(light.color);
                var _light = new THREE.PointLight(color);
                var element = new Float32Array(light.init_orientation);
                var initMatrix = new THREE.Matrix4();
                initMatrix.elements = element;
                _light.applyMatrix(initMatrix);
                break;
            //TODO add other light cases..
            case "SomeOtherLight":
                break;
        } 
        console.log("in Light");
        _light.name = light.simulation_id;
        console.log("in Light");
        _light["object-info"] = light;
        console.log("in Light");
        self._scene.add(_light);

    },

    runAnimation: function(){
        var self = this;
        jQuery("#playAnimation").css("display","none");
        jQuery("#stopAnimation").css("display","block");
        var currentTime = 0;
        var timeDelta = self.model.timeDelta;
        
        self.animationID = window.setInterval(function(){ 
            // setAnimationTime sets slider and scene
            // to that particular time.
            self.setAnimationTime(currentTime);
            currentTime+=timeDelta;

            
        }, timeDelta*1000);

    },

    setAnimationTime: function(currentTime){
        var self = this;
        // Set the slider to the current animation time..
        if(currentTime>=self._finalTime) {
            self.stopAnimation();
        }    
        var percent = currentTime/self._finalTime*100;
        // Now animate objects in scene too..
        var i = self._timeArray.indexOf(currentTime);
        var _children = self._scene.children;

        for(var i=0;i<_children.length;i++){
            var id = _children[i].name;
            if(self.simData[id] != undefined){
                var element = new Float32Array(self.simData[id][i]);
                var orientationMatrix = new THREE.Matrix4();
                orientationMatrix.elements = element;
                //_children[i].matrix.identity()
                _children[i].applyMatrix(orientationMatrix);
            }

        }
        jQuery("#timeSlider").slider("setValue",percent);
        jQuery("#time").html(" " + Math.round(currentTime*100)/100 + "s");
        
    },

    stopAnimation: function(){
        var self = this;
        console.log("INFO: Stopping Animation");
        window.clearInterval(self.animationID);
        self.setAnimationTime(0)
        jQuery("#stopAnimation").css("display","none");
        jQuery("#playAnimation").css("display","block");

    },

    applySceneInfo: function(){
        // TODO.. save data from the objectDialog
        // into the _meshes, as well into codeMirror's JSON..
        // 

        

        
        alert("HEYA");    

    },

    _removeAll: function(){
        var self = this;
        var _children = self._scene.children;

        for(var i=_children.length-1;i>=0;i--) { 
            if(_children[i].name){
                self._scene.remove(_children[i]);
            }
        };



    }




});


