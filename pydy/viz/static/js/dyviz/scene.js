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
        /*
        self._addAxes();
        self._addDefaultLightsandCameras();
        self._addObjects();
        self._render();
        */

	},

	_createRenderer: function(){
		/**
		  * Creates a webGL Renderer
		  * with a default background color.
		  *
		**/ 
		var self = this;
		self.renderer = new THREE.WebGLRenderer();
	    self.renderer.setSize(800, 800);
	    var backgroundColor = new THREE.Color(161192855); // WhiteSmoke
	    self.renderer.setClearColor(backgroundColor);	
	    var container = $('#renderer');
	    container.append(this.renderer.domElement);	

	    // new Scene..


	},

	loadSimulation: function(){

		/** 
		  * This method loads the simulation data 
		  * from the simulation JSON file. The data is
		  * saved in the form of 4x4 matrices mapped to 
		  * the simulation object id, at a particular time.
		  * 
		**/ 
		var self = this;

        var path = self.getBasePath() + self.model.simulationData;
        
        new Ajax.Request(path, {
            method:'get',
            onSuccess: function(transport) {
            	// Got file here.. load this on Canvas!
            	self.simData = $.parseJSON(transport.responseText);

            },
            onFailure: function() { alert('Simulation File not loaded!'); },
            on404: function(){ alert("Simulation File Not Found! Error:404"); }

        }); 


        // TODO: form data into a better Data Structure
        // and use it in the animation

        

	},


});
})(jQuery);

