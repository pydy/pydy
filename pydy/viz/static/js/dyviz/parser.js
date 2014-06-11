
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
		vara filePath = self._sceneFilePath;
		if(filePath.split(".").pop() !== "json"){
			console.log("ALERT: File should be a valid JSON file!");
			alert("File should be a valid JSON file!");
			return;
		}
		
        
        var self = this;
        self._model = {};
        alert(filePath);

	},

	loadSimulation: function(){

		/** 
		  * This method loads the simulation data 
		  * from the simulation JSON file. The data is
		  * saved in the form of 4x4 matrices mapped to 
		  * the simulation object id, at a particular time.
		  * 
		**/ 

        

	},


});
})(jQuery);

