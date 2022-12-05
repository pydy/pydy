
DynamicsVisualizer.Scene = Object.extend(DynamicsVisualizer, {

    create: function(){
        /**
          * This method creates the scene from the self.model
          * and renders it onto the canvas.
        **/
        var self = this;
        self._createRenderer();
        self._createEmptyScene();
        self._addDefaultCamera();
        self._addDefaultLight();
        self._addAxes();
        self._addTrackBallControls();
        self.WindowResize(self.webgl_renderer, self.currentCamera, this);
        self.animationPaused = false;
    },

    _createRenderer: function(){
        /**
          * Creates a webGL Renderer
          * with a default background color.
        **/
        var self = this;
        self.webgl_renderer = new THREE.WebGLRenderer();
        self._updateWidth();
        self._updateHeight();
        self.webgl_renderer.setSize(self.width, self.height);
        var backgroundColor = new THREE.Color(161192855); // WhiteSmoke
        self.webgl_renderer.setClearColor(backgroundColor);
        var container = jQuery('#renderer');
        container.append(self.webgl_renderer.domElement);
    },

    _createEmptyScene: function(){
        /**
          * Creates a THREE Scene
          *
        **/
        var self = this;
        // new Scene..
        self._scene  = new THREE.Scene();
    },

    _addDefaultCamera: function(){
        /**
          * This method adds
          * a Perspective camera to the
          * initial visualization
        **/
        var self = this;

        self.primaryCamera = new THREE.PerspectiveCamera();
        self.primaryCamera.position.set(0,0,100);
        self._updateWidth();
        self._updateHeight();
        self.primaryCamera.aspect = self.width / self.height;
        self._scene.add(self.primaryCamera);
        self.currentCamera = self.primaryCamera;
        self.currentCamera.updateProjectionMatrix();
    },

    _addDefaultLight: function(){
        /**
          * This method adds a default light
          * initial visualization
        **/
        var self = this;
        var light = new THREE.PointLight(0xffffff);
        light.position.set(10,10,-10);
        self._scene.add(light);
    },


    _addAxes: function(){
        /**
          * Adds a default system of axes
          * to the initial visualization.
        **/
        var self = this;

        var axes = new THREE.AxisHelper(100);
        self._scene.add(axes);

    },

    _addTrackBallControls: function(){
        /**
          * Adds Mouse controls
          * to the initial visualization
          * using TrackballControls Library.
        **/
        var self = this;
        self.primaryControls = new THREE.TrackballControls(self.currentCamera,
                                            self.webgl_renderer.domElement);

    },

    _updateWidth: function(){
        var self = this;
        // Setting minimum width to be 800px
        if(jQuery(window).width() > 800) {
            self.width = jQuery(window).width() * 0.69;
        } else{
            self.width = 800 * 0.665;
        }
    },

    _updateHeight: function(){
        var self = this;
        self.height = jQuery(window).height() * 0.83;
    },

    _map_points_to_curve: function(points){
        /**
          * Maps points to a THREE.Curve object.
        **/
        var self = this;
        var vector_array = [];

        // NOTE: Due to some issue in THREE.Vector3, we need to parseFloat
        // coordinates before passing it to Vector3
        for(var point in points){
            vector_array.push(
                new THREE.Vector3(
                    parseFloat(point[0]), parseFloat(point[1]),
                    parseFloat(point[2])
                )
            );
        }

        curve = new THREE.SplineCurve3(vector_array);
        return curve;
    },

    WindowResize: function(renderer, camera, self){
        /**
          * Adds window resize event listener
          * to renderer and camera and updates
          * them accordingly. This is a modified
          * version of THREEX.WindowResize.js
          * LICENCE: The MIT License (MIT)
          *
          * Copyright (c) 2013 Jerome Etienne
          *
          * Permission is hereby granted, free of charge, to any person obtaining a copy of
          * this software and associated documentation files (the "Software"), to deal in
          * the Software without restriction, including without limitation the rights to
          * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
          * the Software, and to permit persons to whom the Software is furnished to do so,
          * subject to the following conditions:
          *
          * The above copyright notice and this permission notice shall be included in all
          * copies or substantial portions of the Software.
        **/
        var callback  = function(){
          self._updateWidth();
          self._updateHeight();
          // notify the renderer of the size change
          renderer.setSize( self.width, self.height );
          // update the camera
          camera.aspect = self.width / self.height;
          camera.updateProjectionMatrix();
        };
        // bind the resize event
        window.addEventListener('resize', callback, false);
        // return .stop() the function to stop watching window resize
        return {
          /**
            * Stop watching window resize
          **/
          stop : function(){
            window.removeEventListener('resize', callback);
          }
        };
    },

    _resetControls: function(){
        /**
          * Resets the scene camera to
          * the initial values(zoom, displacement etc.)
        **/
        var self = this;
        self.primaryControls.reset();
    },

    addObjects: function(){
        /**
          * Adds the geometries
          * loaded from the JSON file
          * onto the scene. The file is
          * saved as an object in self.model
          * and then rendered to canvas with this
          * function.
        **/
        var self = this;

        self._removeAll(); // Removes old objects first!

        var objects = self.model.objects;
        for(var i in objects) self._addIndividualObject(objects[i]);

    },

    addCameras: function(){
        /**
          * Adds the cameras
          * loaded from the JSON file
          * onto the scene.
        **/
        var self = this;
        var cameras = this.model.cameras;
        for(var i in cameras)  self._addIndividualCamera(cameras[i]);

    },

    addLights: function(){
        /**
          * Adds the Lights
          * loaded from the JSON file
          * onto the scene.
        **/
        var self = this;
        var lights = this.model.lights;
        for(var i in lights)  self._addIndividualLight(lights[i]);
    },

    _addIndividualObject: function(object){
        /**
          * Adds a single geometry object
          * which is taken as an argument
          * to this function.
        **/
        var self = this;
        var type = object.type;

        switch(type) {

            case "Cube":
                var geometry = new THREE.CubeGeometry(
                                  object.length,
                                  object.length,
                                  object.length,
                                  50, 50, 50);
                break;

            case "Box":
                var geometry = new THREE.CubeGeometry(
                                  object.width,
                                  object.height,
                                  object.depth,
                                  50, 50, 50);
                break;
                
            case "Sphere":
                var geometry = new THREE.SphereGeometry(
                                               object.radius, 100);
                break;

            case "Cylinder":
                var geometry = new THREE.CylinderGeometry(object.radius,
                                                          object.radius,
                                                          object.length,100);

                break;

            case "Cone":
                var geometry = new THREE.CylinderGeometry(
                                          object.radius,
                                          object.radius/100,
                                          object.length,
                                          50,50, openEnded=true);
                break;

            case "Circle":
                var geometry = new THREE.CylinderGeometry(object.radius,
                                                          object.radius,
                                                          0,100);
                break;

            case "Plane":
                var geometry = new THREE.PlaneGeometry(
                                          object.length,
                                          object.width,
                                          100);
                break;

            case "Tetrahedron":
                var geometry = new THREE.TetrahedronGeometry(
                                          object.radius);
                break;

            case "Octahedron":
                var geometry = new THREE.OctahedronGeometry(
                                          object.radius);
                break;

            case "Icosahedron":
                var geometry = new THREE.IcosahedronGeometry(
                                          object.radius);
                break;

            case "Torus":
                var geometry = new THREE.TorusGeometry(
                                          object.radius,
                                          object.tube_radius,100
                                          );
                break;

            case "TorusKnot":
                var geometry = new THREE.TorusKnotGeometry(
                                          object.radius,
                                          object.tube_radius,100
                                          );
                break;

            case "Tube":
                var curve = self._map_points_to_curve(object.points);
                var geometry = new THREE.TubeGeometry(curve, 64,
                                                      object.radius, 8, false);
                break;

        }

        var material = self.Materials.getMaterial(object.material);
        material.color = new THREE.Color(object.color);
        var mesh = new THREE.Mesh(geometry, material);
        if(type == 'Plane'){
            mesh.material.side = THREE.DoubleSide;
        }
        var element = new Float32Array(object.init_orientation);
        var initMatrix = new THREE.Matrix4();
        initMatrix.elements = element;
        mesh.matrix.identity();
        mesh.applyMatrix(initMatrix);
        mesh["object-info"] = object;
        mesh.name = object.simulation_id;
        self._scene.add(mesh);
    },

    _addIndividualCamera: function(camera){
        /**
          * Adds a single camera object
          * which is taken as an argument
          * to this function.
        **/
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

            case "OrthoGraphicCamera":
                var _camera = new THREE.OrthographicCamera(
                                                -320, 320,
                                                240, -240,
                                                camera.near, camera.far );
                var _element = new Float32Array(camera.init_orientation);
                var initMatrix = new THREE.Matrix4();
                initMatrix.elements = _element;
                _camera.applyMatrix(initMatrix);

        }
        _camera.name = camera.simulation_id;
        _camera["object-info"] = camera;
        self._updateWidth();
        self._updateHeight();
        _camera.aspect = self.width / self.height;
        self._scene.add(_camera);
        self.currentCamera = _camera;
        self.currentCamera.updateProjectionMatrix();
        self._addTrackBallControls();
    },

    _addIndividualLight: function(light){
        /**
          * Adds a single light object
          * which is taken as an argument
          * to this function.
        **/

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
            //TODO some other Light implementations
        }

        _light.name = light.simulation_id;
        _light["object-info"] = light;
        self._scene.add(_light);
    },

    runAnimation: function(){
        /**
          * This function iterates over the
          * the simulation data to render them
          * on the canvas.
        **/
        var self = this;
        // toggle buttons..
        jQuery("#play-animation").prop('disabled', true);
        jQuery("#pause-animation").prop('disabled', false);
        jQuery("#stop-animation").prop('disabled', false);

        var startTime = self.model.startTime;
        if(!self.animationPaused){
          self.currentTime = startTime;
        }

        self.animationPaused = false;
        var timeDelta = self.model.timeDelta;

        self.animationID = window.setInterval(function(){
                self.setAnimationTime(self.currentTime);
                self.currentTime += timeDelta;
                if(self.currentTime >= self._finalTime){
                  self.currentTime = startTime;
                  if(!jQuery("#play-looped").is(":checked")){
                    self.stopAnimation();
                  }
                }
            },
        timeDelta*1000);
    },

    setAnimationTime: function(currentTime){
        /**
          * Takes a time value as the argument
          * and renders the simulation data
          * corresponding to that time value.
        **/
        var self = this;
        var t0 = self.model.startTime;
        var percent = (100*(currentTime - t0)/(self._finalTime - t0)).toFixed(3);

        var time_index = self._timeArray.indexOf(currentTime);
        var _children = self._scene.children;
        for(var i=0;i<_children.length;i++){
          if(!(_children[i] instanceof (THREE.OrthoGraphicCamera || THREE.PerspectiveCamera))){
            var id = _children[i].name;
            if(self.simData[id] != undefined){
                var element = new Float32Array(self.simData[id][time_index]);
                var orientationMatrix = new THREE.Matrix4();
                orientationMatrix.elements = element;
                _children[i].matrix.identity()
                _children[i].applyMatrix(orientationMatrix);
            }
          }
        }
        jQuery("#time-slider").slider("setValue", percent);
        jQuery("#time").html(" " + currentTime.toFixed(3) + " s");

    },

    pauseAnimation: function(){
       /**
         * Pauses the animation at the
         * current frame.
       **/
       var self = this;
       console.log("[PyDy INFO]: Pausing Animation");
       jQuery("#play-animation").prop('disabled', false);
       jQuery("#pause-animation").prop('disabled', true);
       jQuery("#stop-animation").prop('disabled', false);
       window.clearInterval(self.animationID);
       self.animationPaused = true;

    },
    stopAnimation: function(){
        /**
          * Stops the animation, and
          * sets the current time value to initial.
        **/
        var self = this;
        console.log("[PyDy INFO]: Stopping Animation");
        if(!self.animationPaused){
          window.clearInterval(self.animationID);
        }
        self.currentTime = self.model.startTime;
        self.setAnimationTime(self.currentTime);
        jQuery("#play-animation").prop('disabled', false);
        jQuery("#pause-animation").prop('disabled', true);
        jQuery("#stop-animation").prop('disabled', true);

    },

    _removeAll: function(){
        /**
          * Removes all the geometry elements
          * added to the scene from the loaded scene
          * JSON file. Keeps the default elements, i.e.
          * default axis, camera and light.
        **/
        var self = this;
        var _children = self._scene.children;

        for(var i=_children.length-1;i>=0;i--) {
            if(_children[i].name){
                self._scene.remove(_children[i]);
            }
        };
    },

    _blink: function(id){
        /**
          * Blinks the geometry element.
          * takes the element simulation_id as the
          * argument and blinks it until some event is
          * triggered(UI button press)
        **/
        var self = this;
        self._blinker = self._scene.getObjectByName(id);
        console.log("BLinker: " + self._blinker.name)
        self._blinker.visible = false;
        self.blinkId = window.setInterval(function(){
            if(self._blinker.visible == false){
                self._blinker.visible = true;
            } else{
                self._blinker.visible = false;
            }
        }, 500);
    }
});
