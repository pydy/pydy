
//Modifying AnimFrame for smart animations for different browsers
window.requestAnimFrame = (function(){
      return  window.requestAnimationFrame       || 
              window.webkitRequestAnimationFrame || 
              window.mozRequestAnimationFrame    || 
              window.oRequestAnimationFrame      || 
              window.msRequestAnimationFrame     || 
              function(/* function */ callback, /* DOMElement */ element){
                window.setTimeout(callback, 1000 / 60);
              };
    })();	
    
    
    




// fetching coorinates from Python output
var coordinates = #(coordinates_time);

var link_length = coordinates[0][0];
var time = coordinates[1][0];
var i=0;
//Setting other varaibles
var controls,scene,camera,renderer,reset; 
var $reset_button;
var $canvas
var sphere_geometry;

var sphere, link;
var animation_request_id;


//initial Position, dimensions of bob
var radius = 0.1, segments = 100, rings = 100;

init_canvas();
init();
visualize();






	
function init_canvas(){
	
	//This function initiates a basic canvas, with trackball controls, 
	//all drawing work occurs in init() function.
	
	// first of all , a renderer ...
	renderer = new THREE.WebGLRenderer();
	
	//show the IPython handle container
	container.show();
	
	// create a canvas div
	$canvas = $("<div/>").attr("id", "#canvas");
	//giving background color
	
	$canvas.attr("style","background-color:rgb(104,104,104)");
	
	//For this particular, giving a top left, position time index
	
	$canvas.append("<p>click left and move mouse to rotate camera, hit reset button to reset camera</p>");
	$canvas.append('<div id=\"pos_index\" style=\"margin-left:10px;margin-top:20px;\">T= ' + time + '<br /> Pos..='+(-link_length)+'</div>');
	
	// Adding our canvas to IPython UI
	
	
	// Now lets add a scene ..
	
	// set the scene size
	 var WIDTH = 400,
	    HEIGHT = 300;
	    
	scene = new THREE.Scene();
	
	
	//Add a camera to the scene..
	
	// set some camera attributes
	var VIEW_ANGLE = 45,
	    ASPECT = WIDTH / HEIGHT,
	    NEAR = 0.1,
	    FAR = 1000;
	    
	camera = new THREE.PerspectiveCamera(  VIEW_ANGLE,
	                                ASPECT,
	                                NEAR,
	                                FAR  );
	    
	        
	// the camera starts at 0,0,0 so pull it back
	camera.position.z = 10;
	
	
	// Add trackball controls
	
	controls = new THREE.TrackballControls(camera, renderer.domElement);
	
	
	reset = function(){ controls.reset();}
	
	
	scene.add(camera);
	
	

    $reset_button = $('<button/>').attr('style','margin-left:40px;').click(reset);
    $reset_button.append('Reset Camera');                             
    $canvas.append($reset_button)
	
	$canvas.append(renderer.domElement);
	element.append($canvas);
	// start the renderer
	renderer.setSize(WIDTH, HEIGHT);
	
	
	// Add  axes ...
	
	var axesMaterial = new THREE.MeshLambertMaterial(
	{
	    color: 0xFFFFFF
	    
	});
	var x_axis = new THREE.Mesh(
	   new THREE.CubeGeometry(WIDTH, 0.03, 0.03),
	   axesMaterial);
	   
	scene.add(x_axis);
	
	var y_axis = new THREE.Mesh(
	   new THREE.CubeGeometry(0.03, HEIGHT, 0.03),
	   axesMaterial);
	   
	scene.add(y_axis);
	
	var z_axis = new THREE.Mesh(
	   new THREE.CubeGeometry(0.03, 0.03, 10),
	   axesMaterial);
	   
	scene.add(z_axis);
	
	// create a point light
	var pointLight = new THREE.PointLight( 0xFFFFFF );

	// set its position
	pointLight.position.x = 100;
	pointLight.position.y = 100;
	pointLight.position.z = 100;

	// add to the scene
	scene.add(pointLight);
	
	renderer.render(scene, camera);    
	
	
}
	




function init(){
	
    // create the sphere's material
	var sphereMaterial = new THREE.MeshLambertMaterial(
	{
	    color: 0xCC0000
	});

	// create a new mesh with sphere geometry
	
	sphere_geometry = new THREE.SphereGeometry(radius, segments, rings)
	sphere = new THREE.Mesh(sphere_geometry
	   ,
	   sphereMaterial);
    sphere.position.y = -link_length;
    sphere.position.needsUpdate = true;
    sphere.geometry.dynamic = true;
	
	// add the sphere to the scene
	scene.add(sphere);
    	
    // draw!
	renderer.render(scene, camera);    
    
    }



function visualize() {
	
	
	controls.update();
	
	
    
    renderer.render(scene, camera);
	
	requestAnimationFrame(visualize);    
	
   
}

var animate  =function() { 
	
	
	
	link_length = coordinates[0][i];
	time = coordinates[1][i];
	console.log('link_length = '+link_length)
	sphere.position.y = -link_length;
	
	
	console.log(sphere.position.y);
	console.log('radii = '+sphere.geometry.radius)
	$("#pos_index").html('T= ' + time + '<br /> Pos..='+(-link_length));
	
	renderer.render(scene, camera);
	i+=1;
	if(i == coordinates[0].length)
		i=0;
		
	
	animation_request_id = requestAnimationFrame(animate);    
	}	
	
var stop_animate  =function() { 
	
	
	cancelAnimationFrame(animation_request_id);
	 
	}		
	
	
var $animate_button = $('<button/>').attr('style','margin-left:40px;').click(animate);
$animate_button.append('Animate');                             
$canvas.append($animate_button)

var $stop_animate_button = $('<button/>').attr('style','margin-left:40px;').click(stop_animate);
$stop_animate_button.append('Stop Animation');                             
$canvas.append($stop_animate_button)

