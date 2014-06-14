(function($) {

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
            //self._render_with_trackball(); // render with trackball controls
            
            /* should work when load simulation button is pressed!
            self._addObjects();
            */
        
    	},

    	_createRenderer: function(){
    		/**
    		  * Creates a webGL Renderer
    		  * with a default background color.
    		  *
    		**/ 
    		var self = this;
    		self.webgl_renderer = new THREE.WebGLRenderer();
    	    self.webgl_renderer.setSize(800, 600);
    	    var backgroundColor = new THREE.Color(161192855); // WhiteSmoke
    	    self.webgl_renderer.setClearColor(backgroundColor);	
    	    var container = $('#renderer');
    	    container.append(self.webgl_renderer.domElement);	
            
            // new Scene..
            self._scene  = new THREE.Scene();

        },    
            
    	_addDefaultLightsandCameras: function(){
    		var self = this;
            self._cameras = {};

    		var camera = new THREE.PerspectiveCamera();
            camera.position.x = 00 
            camera.position.y = 00 
            camera.position.z = 100
            self._cameras["init_camera"] = camera;
            console.log(self._cameras.init_camera);
            self._scene.add(self._cameras.init_camera);

            self._lights = {};
            var light = new THREE.PointLight(0xffffff);
            light.position.set(10,10,-10);
            self._lights["init_light"] = light;
            console.log(self._lights.init_light)
            self._scene.add(self._lights.init_light);

    	},

    	_addAxes: function(){
    		var self = this;
    		var axes = new THREE.AxisHelper(100);
            axes.name = "Axes";
            self._scene.add(axes);
            
    	},

        _addTrackBallControls: function(){
        	
        	this.primaryControls = new THREE.TrackballControls(this._cameras.init_camera,
                                                this.webgl_renderer.domElement);
        
        },

        _resetControls: function(){
        	this.primaryControls.reset();
        },

        /* Doesnt work.. not sure why!! 
        _render_with_trackball: function(){
        	console.log("In this function");
        	var self = this;
        	var renderer = self.webgl_renderer;
        	requestAnimationFrame(this._render_with_trackball);
        	
        	var self = this;
        	this.webgl_renderer.render(self._scene, self._camera);
        	self.primaryControls.update();
        },
        */

        _addObjects: function(){
            var self = this;
            var objects = self.model.objects;
            self._frames = {};
            self._cameras = {};
            self._lights = {};
            for(var i=0;i<objects.length; i++) self._scene.add(objects[i]);
        },




    });
})(jQuery);


