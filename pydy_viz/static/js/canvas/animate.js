
Canvas.prototype.startAnimation = function(){
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
      {
          var _element = JSONObj.frames[key].simulation_matrix[counter];
          var matrix = new THREE.Matrix4();
          matrix.elements = _element;

          Canvas.prototype.scene.getObjectByName("frames").children[key].matrix.identity();
          Canvas.prototype.scene.getObjectByName("frames").children[key].applyMatrix(matrix)
          
       };


Canvas.prototype.pauseAnimation = function(){
    cancelAnimationFrame(Canvas.prototype.animationID);
    Canvas.prototype.animationID = undefined;
    $("#startAnimation").click(Canvas.prototype.startAnimation);  
        
};
Canvas.prototype.stopAnimation = function(){
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

