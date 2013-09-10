
Canvas.prototype.switchCamera = function(){

    Canvas.prototype.cameraCounter++;
    if(Canvas.prototype.cameraCounter >=Canvas.prototype.cameras.children.length)
    {
        Canvas.prototype.cameraCounter = 0;
    }      
    var _camera = Canvas.prototype.cameras.children[Canvas.prototype.cameraCounter];
    renderer.render(scene, _camera);
    console.log("Switched to camera:" + Canvas.prototype.cameraCounter);    
   
};


Canvas.prototype.goToFrame = function(frame){ var i = $("#frame-number").val(); alert(i); }


