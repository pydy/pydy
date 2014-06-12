
var DynamicsVisualizer = new DynamicsVisualizer();

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
		self._model = {};


		if(self.getFileExtenstion() !== "json"){
			console.log("ALERT: File should be a valid JSON file!");
			alert("File should be a valid JSON file!");
			return;
		}

		new Ajax.Request(filePath, {
            method:'get',
            onSuccess: function(transport) {
            	// Got file here.. load this on Canvas!
            	self.model = $.parseJSON(transport.responseText);

            },
            onFailure: function() { alert('Scene File not loaded!'); },
            on404: function(){ alert("Scene File Not Found! Error:404"); }

        });
		
		// loads scene onto canvas!
		//self.Scene.load();
		

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

