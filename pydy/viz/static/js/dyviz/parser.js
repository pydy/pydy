
(function($) {

	DynamicsVisualizer.Parser = Object.extend(DynamicsVisualizer, {
		
		loadScene: function(){

			/** 
			  * This method calls an ajax request on the 
			  * JSON file and reads the scene info from 
			  * the JSON file. 
			  * 
			  * 
			**/ 
			var self = this;
			var filePath = self.sceneFilePath;
			
			if(self.getFileExtenstion() !== "json"){
				console.log("ALERT: File should be a valid JSON file!");
				alert("File should be a valid JSON file!");
				return;
			}

			new Ajax.Request(filePath, {
	            method:'get',
	            onSuccess: function(transport) {
	            	// Got file here.. load this on Canvas!
	            	// This call is not complete, before the Scene.
	            	self.model = $.parseJSON(transport.responseText);

	            },
	            onComplete: function(){
	            	console.log("request completed, adding Objects to scene");
	            	self.Scene.addObjects();
	            	self.Scene.addCameras();
	                self.Scene.addLights();
	                self.Parser.loadSimulation();
	                // activate run Animation button!
	                $("#playAnimation").removeClass("disabled");
	                
	                

	 
	            },
	            onFailure: function() { alert('Scene File not loaded!'); },
	            on404: function(){ alert("Scene File Not Found! Error:404"); }

	        });
			

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
	        console.log(path);
	        new Ajax.Request(path, {
	            method:'get',
	            onSuccess: function(transport) {
	            	self.simData = $.parseJSON(transport.responseText);
	            	


	            },
	            onComplete: function(){
	            	// build timeArray and dataArray from simdata...
	            	
	            	self.createTimeArray();
	            	//self.createDataArray();
	            },

	            onFailure: function() { alert('Simulation File not loaded!'); },
	            on404: function(){ alert("Simulation File Not Found! Error:404"); }

	        }); 

		},

		createTimeArray: function(){
			var self = this;
			var _NtimeSteps = self.model.timeSteps;
			var timeDelta = self.model.timeDelta;
			var time=0;
			self._timeArray = [];

			for(var i=0;i<_NtimeSteps; i++){
				self._timeArray.push(time);
				time+=timeDelta;
			}
			self._finalTime = self._timeArray.slice(-1)[0];
			console.log("Created Time Array");
		}


	});
})(jQuery);

