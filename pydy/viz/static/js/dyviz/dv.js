// Dynamics Visualizer main class.

var DynamicsVisualizer = {};

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
			jQuery('head').append('<link rel="stylesheet" type="text/css" href="css/main.css">');
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
		jQuery("#simulation-load").click(function(){
			// Fetch FilePath from input..
			// load it into the canvas

            self.sceneFilePath = jQuery("#json-input").val();
            console.log("INFO: found scene JSON file:" + self.sceneFilePath);
			self.Parser.loadScene();

		});

		jQuery("#json-save").click(function(){
			// Activate CodeMirror... 
			//
		});

		jQuery("#timeSlider").slider({min:0,max:100,step:1, handle:"square", value:0});

		jQuery("#resetControls").click(function(){
			// Activate CodeMirror... 
			//
		});

		jQuery("#playAnimation").click(function(){
			self.Scene.runAnimation();
			
		});
		jQuery("#stopAnimation").click(function(){
			self.Scene.stopAnimation();
			
		});
		
		console.log("INFO: Activated UI controls");


	},

	loadUIElements: function(){
		var self = this;
        console.log("Here1");
        jQuery("#playAnimation").removeClass("disabled");

        var objs = self.model.objects;
        console.log("Here:" + objs);
        for(var obj=0;obj<objs.length;obj++){

            var toAppend = '<li><a id="'+ objs[obj].simulation_id + 
                           '" href="#">' + objs[obj].name + '</a></li>';

            jQuery("#objectDropdown").append(toAppend);

            // adding click functions to all dropdown objs.

            jQuery("#" + objs[obj].simulation_id).click(function(){
            	self._openDialog(jQuery(this).attr("id"));
          
            });

        }
    },   

    _openDialog: function(obj){
    	var self = this;
    	console.log("In Open dialog:")
    	var toLoad = self._meshes[obj]["object-info"];

    	var mainDiv = jQuery('<div/>',{id: "object-"+ toLoad.simulation_id, style: 'display:none;'}); 
    	
    	// for name..
        var div1 = jQuery('<div />',{class: 'input-group'});
    	div1.append('<span class="input-group-addon">Name</span>');
    	div1.append(jQuery('<input />',{ type:'text', id: "_name", class: 'form-control', value: toLoad.name}));

        // for color..
    	div1.append('<span class="input-group-addon">Color</span>');
    	div1.append(jQuery('<input />',{ type:'text', id: "_color", class: 'form-control', value: toLoad.color}));
        
        // for material.. a dropdown..

        var div_geom = jQuery('<select />',{class: 'form-control'});
        div_geom.append('<option value="' + toLoad.material + '">' + toLoad.material + '</option>');
        for(var i in self.Materials){
        	div_geom.append('<option value="' + i +  '">' + i + '</option>');
        }
        
        //rest geom. params depend on shape..
        var div2 = jQuery('<div />',{class: 'input-group'});

        switch(toLoad.type){
        	// TODO all objects..
        	case "Sphere":
        	    div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Radius'));
    	        div2.append(jQuery('<input />',{ type:'text', id: "_radius", class: 'form-control', value: toLoad.radius}));
        	    break;

        	case "Cylinder":
        	    div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Radius'));
    	        div2.append(jQuery('<input />',{ type:'text', id: "_radius", class: 'form-control', value: toLoad.radius}));

    	        div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Length'));
    	        div2.append(jQuery('<input />',{ type:'text', id: "_length", class: 'form-control', value: toLoad.length}));
        	    break;

        }





        mainDiv.append(div1);
        mainDiv.append('<hr/><span class="input-group-addon">Material</span>');
        mainDiv.append(div_geom);
        mainDiv.append("<hr />");
        mainDiv.append(div2);
        // finally a button..
        mainDiv.append('<hr /><button id="scene-info-apply" class="btn btn-primary btn-small">Apply</button>');
        jQuery("#objectDialog").html(mainDiv);

        // show after whole div is populated..
        mainDiv.fadeIn("slow");

        // finally activate button..
        jQuery("#scene-info-apply").click(function(){
			self.Scene.applySceneInfo();
			
		});

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

var DynamicsVisualizer = new DynamicsVisualizer();


