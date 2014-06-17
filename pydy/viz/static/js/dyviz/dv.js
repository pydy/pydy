// Dynamics Visualizer main class.

var DynamicsVisualizer = {};

(function($) {

DynamicsVisualizer = Class.create({
	/** 
	  * DV is the main class for Dynamics Visualizer.
	  * 
	**/ 

	init: function(){    
		/**
		  * Initializes the Class with 
		  * relevant objects.
		  *
		**/
         
		var self = this;
		console.log("INFO: initializing Visualizer");
		

		this._model = {}; 
        
		if(!self.isWebGLCompatible()){
			console.log("ALERT: Incompatible browser!");
			alert("The browser you are using is not compatible! " + 
				"Please use a latest version of Chrome or Firefox");

		    return false;
		}
		// Load different css for Browser and IPython notebook...
		/*try{
			jQuery('head').append('<link rel="stylesheet" type="text/css" href="css/ipython_main.css">');
		}
		catch(err){
			alert("IPython notebook not triggered!")
			console.log("IPython notebook not triggered, loading main.css");
			$('head').append('<link rel="stylesheet" type="text/css" href="css/main.css">');
		}*/

    },


	isWebGLCompatible: function(){ 
		/**
		  * Checks whether the browser used is
		  * compatible for handling webGL based
		  * animations. Raises an alert if not!
		  * 
		  * Requires external script: Modernizr.js
		  *
		**/
        
        if (!Modernizr.canvas || !Modernizr.webgl) return false;
        else return true;


	   },

	activateUIControls: function(){
		/**
		  * This method adds functions to the UI buttons
		  * It should be **strictly** called after the 
		  * other DynamicsVisualizer sub-modules are loaded
		  * in the browser, else certain functionality will 
		  * be(not might be!) hindered.
		**/


        var self = this;
		$("#simulation-load").click(function(){
			// Fetch FilePath from input..
			// load it into the canvas

            self.sceneFilePath = $("#json-input").val();
            console.log("INFO: found scene JSON file:" + self.sceneFilePath);
			self.Parser.loadScene();

		});

		$("#json-save").click(function(){
			// Activate CodeMirror... 
			//
		});

		$("#timeSlider").slider({min:0,max:100,step:1, handle:"square", value:0});

		$("#resetControls").click(function(){
			// Activate CodeMirror... 
			//
		});

		$("#playAnimation").click(function(){
			self.Scene.runAnimation();
			
		});
		$("#stopAnimation").click(function(){
			self.Scene.stopAnimation();
			
		});
		$("#distractionFreeMode").click(function(){
			self.Utils.distractionFreeMode();
			
		});


		console.log("INFO: Activated UI controls");


	},

    getBasePath: function(){
    	var self = this;
    	var slashes_fixed = self.sceneFilePath.replace(/\\/g, "/");
        return slashes_fixed.split("/").slice(0,-1).join("/") + "/";
    },

    getFileExtenstion: function(){
    	var self = this;
    	return self.sceneFilePath.split(".").slice(-1)[0].toLowerCase();

    }

    });

    
})(jQuery);


var DynamicsVisualizer = new DynamicsVisualizer();


