
DynamicsVisualizer.Parser = Object.extend(DynamicsVisualizer, {

    loadScene: function(){

        /**
          * This method calls an ajax request on the
          * JSON file and reads the scene info from
          * the JSON file, and saves it as an object
          * at self.model.
        **/
        var self = this;
        var filePath = self.sceneFilePath;
        if(self.getFileExtenstion() !== "json"){
            alert("[PyDy ALERT]: File should be a valid JSON file!");
            return;
        }

        new Ajax.Request(filePath, {
            method:'get',
            onSuccess: function(transport) {
                // Got file here.. load this on Canvas!
                // This call is not complete, before the Scene.
                self.model = jQuery.parseJSON(transport.responseText);

            },
            onComplete: function(){
                console.log("[PyDy INFO]: Ajax request completed, adding Objects to scene");
                self.Scene.addObjects();
                self.Scene.addCameras();
                self.Scene.addLights();
                self.Parser.loadSimulation();
                // Load UI elements relevant for animation stuff!
                self.loadUIElements();
                console.log("Done with loadScene");
            },
            onFailure: function() {
                alert('[PyDy ALERT]: Scene file not loaded!');
                console.log('[PyDy ALERT]: Scene file not loaded!');
            },
            on404: function(){
                alert("[PyDy ALERT]: Scene file not Found! error:404");
                console.log("[PyDy ALERT]: Scene file not found! error:404");
            }
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
                self.simData = jQuery.parseJSON(transport.responseText);
            },
            onComplete: function(){
                self.createTimeArray();
                self.Scene.setAnimationTime(self.model.startTime);
            },
            onFailure: function() { alert('[PyDy ALERT]: Simulation File not loaded!'); },
            on404: function(){ alert("[PyDy ALERT]: Simulation File Not Found! Error:404"); }
        });
    },

    createTimeArray: function(){
        /**
          * Creates a time array from
          * the information inferred from
          * simulation data.
        **/
        var self = this;
        var timeSteps = self.model.timeSteps;
        var timeDelta = self.model.timeDelta;
        var time = self.model.startTime;
        self._timeArray = [];

        for(var i = 0; i < timeSteps; i++){
            self._timeArray.push(time);
            time += timeDelta;
        }
        self._finalTime = self._timeArray.slice(-1)[0];
        console.log("[PyDy INFO]: Created Time Array");
    }
});
