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

		self._slider = jQuery("#timeSlider").slider({min:0,max:100,step:1, handle:"square", value:0});
        jQuery("#timeSlider").fadeOut();
        self._slider.on('slideStop',function(ev) { 
        	var val = ev.value;
        	var len = self._timeArray.length;
        	var i = 0;
        	var gotValue = false;
        	for( i=0;i<self._timeArray.length && !gotValue;i++){
        		
        		var percent = (self._timeArray[i]/self._timeArray[len-1])*100;
        		if(val <= percent){
        			gotValue = true;
        			break;
        		}
        	//
            }
            self.Scene.setAnimationTime(self._timeArray[i]);

        });
			
			
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
		jQuery("#close-object-dialog").click(function(){
			jQuery("#objectDialog").html(" ");
			jQuery(this).addClass("disabled");

		});

        console.log("INFO: Activated UI controls");


	},

	loadUIElements: function(){
		var self = this;
        
        jQuery("#playAnimation").removeClass("disabled");
        var objs = self.model.objects;
        // clear dropdown of some old objects..
        jQuery("#objectDropdown").find("li").remove()
        
        for(var obj in objs){

            var toAppend = '<li><a id="'+ objs[obj].simulation_id + 
                           '" href="#">' + objs[obj].name + '</a></li>';

            jQuery("#objectDropdown").append(toAppend);

            // adding click functions to all dropdown objs.

            jQuery("#" + objs[obj].simulation_id).click(function(){
            	self.ParamEditor.openDialog(jQuery(this).attr("id"));
          
            });

        }

        // TODO add constant map to sim-params box..

        var constants = self.model.constant_map;
        console.log(constants);
        var div = jQuery("#simulation-params").fadeOut();
        for(var i in constants){
        	alert("Here i:" + i);
        	div.append('<span class="input-group-addon">' + i + '</span>');
        	div.append(jQuery('<input />',{ type:'text', id: i, class: 'form-control', value: constants[i]}));
        }
        div.fadeIn();
        
        

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


