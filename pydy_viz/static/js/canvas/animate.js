
Canvas.prototype.startAnimation = function(){
/**
 * This prototype function kick starts the animation.
 * It iterates over the frames and apply transformation matrices
 * from Simulation Matrix of that frame, iteratively.
 * by default animation is done for a single loop,
 * which can be changed to looped by the check button from the UI.
 */    
    $("#startAnimation").unbind("click");
        
    for(var key in JSONObj.frames){
      Canvas.prototype._animate(key, Canvas.prototype.animationCounter);
      }
      Canvas.prototype.animationCounter++; 
      $("#currentFrame").html("Current Frame: " + Canvas.prototype.animationCounter);
      
      Canvas.prototype.animationID = requestAnimationFrame(Canvas.prototype.startAnimation);
      if(Canvas.prototype.animationCounter >= Canvas.prototype.timeSteps)
          if($("#isLooped").is(':checked')) Canvas.prototype.animationCounter = 0; 
          else { 
          cancelAnimationFrame(Canvas.prototype.animationID);          
          Canvas.prototype.animationCounter = 0;  
          }         

};

Canvas.prototype._animate = function(key, counter)
/**
 * This prototype function is a helper function for 
 * Canvas.prototype.startAnimation
 * It is used to apply transformation matrices to the frames.
 */
      {
          var _element = JSONObj.frames[key].simulation_matrix[counter];
          var matrix = new THREE.Matrix4();
          matrix.elements = _element;

          Canvas.prototype.scene.getObjectByName("frames").children[key].matrix.identity();
          Canvas.prototype.scene.getObjectByName("frames").children[key].applyMatrix(matrix)
          
       };


Canvas.prototype.pauseAnimation = function(){
/**
 * This prototype function pauses the animation, but retains the
 * current animation frame.
 */    
    cancelAnimationFrame(Canvas.prototype.animationID);
    Canvas.prototype.animationID = undefined;
    $("#startAnimation").click(Canvas.prototype.startAnimation);  
        
};
Canvas.prototype.stopAnimation = function(){
/**
 * This prototype function stops the animation, and resets 
 * current animation frame to 0.
 */    
    Canvas.prototype.animationCounter = 0;
    cancelAnimationFrame(Canvas.prototype.animationID);
    Canvas.prototype.animationID = undefined;
    
    for(var key in JSONObj.frames){
          Canvas.prototype._animate(key, 0);
    }
    $("#startAnimation").click(Canvas.prototype.startAnimation);  
                    
};

Canvas.prototype._transform = function(key, i){
    
    var _element = 
            JSONObj.frames[key].simulation_matrix[i];
    this.frames.children[key].matrix.elements = _element;
    };

