
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
    
    
    


var controls,scene,camera,renderer;
// initial Position, dimensions of link
var link_width = 0.05, final_radii = 0.05 , link_length=#(y);

//initial Position, dimensions of bob
var radius = 0.1, segments = 100, rings = 100;


reset();
function reset(){
	init();

	requestAnimationFrame(animate);    
}
	

	
function init(){
	renderer = new THREE.WebGLRenderer();
	
	
	// Use IPython handles to create a div below the current line
    container.show();
    var $container = $("<div/>").attr("id", "#container");
    
    $container.empty();
    $container.attr("style","background-color:rgb(104,104,104)");
    $container.append("<p style=\"margin-left:10px;margin-top:20px;\">T= #(time)<br /> Pos..=(0,#(y),0)</p>");
    
    $container.append(renderer.domElement);
   // $container.append("<button style=\"margin-top:2px;margin-left:2px;\" onClick=\"reset();\">Reset</button>");
    element.append($container);
	
	// set the scene size
	 var WIDTH = 400,
	    HEIGHT = 300;

	// set some camera attributes
	var VIEW_ANGLE = 45,
	    ASPECT = WIDTH / HEIGHT,
	    NEAR = 0.1,
	    FAR = 1000;

	
	// create a camera..
	// and a scene
	
	
	camera = new THREE.PerspectiveCamera(  VIEW_ANGLE,
	                                ASPECT,
	                                NEAR,
	                                FAR  );
	scene = new THREE.Scene();

	// the camera starts at 0,0,0 so pull it back
	camera.position.z = 10;

	// start the renderer
	renderer.setSize(WIDTH, HEIGHT);

	//Create the axes:
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
	
	var y_axis = new THREE.Mesh(
	   new THREE.CubeGeometry(0.03, 0.03, 10),
	   axesMaterial);
	   
	scene.add(y_axis);
	// create the sphere's material
	var sphereMaterial = new THREE.MeshLambertMaterial(
	{
	    color: 0xCC0000
	});

	// create a new mesh with sphere geometry
	
	
	
	
	
	var sphere = new THREE.Mesh(
	   new THREE.SphereGeometry(radius, segments, rings),
	   sphereMaterial);
    sphere.position.y = -link_length;
    
    sphere.geometry.dynamic = true;
    sphere.geometry.verticesNeedUpdate = true;
    sphere.geometry.normalsNeedUpdate = true;
	// add the sphere to the scene
	scene.add(sphere);
	
	//Add a Plane to represent our spring
		
		
		
		
	
	var link = new THREE.Mesh(new THREE.CylinderGeometry(link_width, final_radii, link_length),sphereMaterial);
    
	link.position.y = -link_length/2 ;
	
	link.geometry.dynamic = true;
    link.geometry.verticesNeedUpdate = true;
    link.geometry.normalsNeedUpdate = true;
	
	scene.add(link);

	// and the camera
	scene.add(camera);

	// create a point light
	var pointLight = new THREE.PointLight( 0xFFFFFF );

	// set its position
	pointLight.position.x = 100;
	pointLight.position.y = 100;
	pointLight.position.z = 100;

	// add to the scene
	scene.add(pointLight);

	//Add  trackball controls
	controls = new THREE.TrackballControls(camera, renderer.domElement);
	
	// draw!
	
    
    
	
	
	
	
	
}



function animate() {
	controls.update();    
	renderer.render(scene, camera);    
	requestAnimationFrame(animate);    
}
	
