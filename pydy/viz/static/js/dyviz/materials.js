DynamicsVisualizer.Materials = {};
DynamicsVisualizer.Geometries = ['Cube', 'Cylinder',
                                 'Cone', 'Sphere', 'Circle',
                                 'Plane', 'Tetrahedron',
                                 'Octahedron', 'Icosahedron',
                                 'Torus', 'TorusKnot', 'Tube',]
//                                 'Mesh'];
DynamicsVisualizer.Materials["default"] = new THREE.MeshLambertMaterial();
DynamicsVisualizer.texture_base_path = "static/textures/"

// If IPython,
if(typeof IPython == "undefined"){

    DynamicsVisualizer.texture_base_path = "textures/";

}


function loadMaterials() {

    var checkerBoardTexture = THREE.ImageUtils.loadTexture( DynamicsVisualizer.texture_base_path + "checkerboard.jpg" );
    checkerBoardTexture.wrapS = checkerBoardTexture.wrapT = THREE.RepeatWrapping;
    checkerBoardTexture.anisotropy = 4;

    var metalTexture = THREE.ImageUtils.loadTexture( DynamicsVisualizer.texture_base_path + "metal.jpg" );
    metalTexture.wrapS = metalTexture.wrapT = THREE.RepeatWrapping;
    metalTexture.repeat.set(2,2);
    metalTexture.anisotropy = 4;

    var moonTexture = THREE.ImageUtils.loadTexture( DynamicsVisualizer.texture_base_path + "moon.jpg" );
    moonTexture.wrapS = moonTexture.wrapT = THREE.RepeatWrapping;
    moonTexture.anisotropy = 4;

    var earthTexture = THREE.ImageUtils.loadTexture( DynamicsVisualizer.texture_base_path + "earth.jpg" );
    earthTexture.wrapS = earthTexture.wrapT = THREE.RepeatWrapping;
    earthTexture.anisotropy = 4;

    var grassTexture = THREE.ImageUtils.loadTexture( DynamicsVisualizer.texture_base_path + "grass.jpg" );
    grassTexture.wrapS = grassTexture.wrapT = THREE.RepeatWrapping;
    grassTexture.repeat.set(1,1);
    grassTexture.anisotropy = 4;

    var dirtTexture = THREE.ImageUtils.loadTexture( DynamicsVisualizer.texture_base_path + "dirt.jpg" );
    dirtTexture.wrapS = dirtTexture.wrapT = THREE.RepeatWrapping;
    dirtTexture.repeat.set(1,1);
    dirtTexture.anisotropy = 4;

    var waterTexture = THREE.ImageUtils.loadTexture( DynamicsVisualizer.texture_base_path + "water.jpg" );
    waterTexture.wrapS = waterTexture.wrapT = THREE.RepeatWrapping;
    waterTexture.repeat.set(1,1);
    waterTexture.anisotropy = 4;

    var lavaTexture = THREE.ImageUtils.loadTexture( DynamicsVisualizer.texture_base_path + "lavatile.jpg" );
    lavaTexture.repeat.set( 4, 2 );
    lavaTexture.wrapS = lavaTexture.wrapT = THREE.RepeatWrapping;
    lavaTexture.anisotropy = 4;

    var shininess = 50, specular = 0x333333, shading = THREE.SmoothShading;

    DynamicsVisualizer.Materials["checkerboard"] = new THREE.MeshPhongMaterial(
        {   map: checkerBoardTexture,
            bumpMap: null,
            bumpScale: 0,
            color: 0xffffff,
            ambient: 0x555555,
            specular: 0x222222,
            shininess: 10,
            shading: shading } );


    DynamicsVisualizer.Materials["metal"] = new THREE.MeshPhongMaterial(
        {   map: metalTexture,
            bumpMap: metalTexture,
            bumpScale: 0.1,
            color: 0xffffff,
            ambient: 0x555555,
            specular: 0x222222,
            shininess: 70,
            shading: shading } );

    DynamicsVisualizer.Materials["dirt"] = new THREE.MeshPhongMaterial(
        {   map: dirtTexture,
            bumpMap: dirtTexture,
            bumpScale: 0.05,
            color: 0xffffff,
            ambient: 0x555555,
            specular: 0x222222,
            shininess: 10,
            shading: shading } );

    DynamicsVisualizer.Materials["foil"] = new THREE.MeshPhongMaterial(
        {   map: waterTexture,
            bumpMap: waterTexture,
            bumpScale: 0.3,
            color: 0xccccaa,
            ambient: 0x555544,
            specular: 0x777777,
            shininess: 100,
            shading: shading } );

    DynamicsVisualizer.Materials["water"] = new THREE.MeshPhongMaterial(
        {   map: waterTexture,
            bumpMap: waterTexture,
            bumpScale: 0.05,
            color: 0x3333aa,
            ambient: 0x335577,
            specular: 0x555555,
            shininess: shininess,
            shading: shading } );


    DynamicsVisualizer.Materials["grass"] = new THREE.MeshPhongMaterial(
        {   map: grassTexture,
            bumpMap: grassTexture,
            bumpScale: 0.05,
            color: 0xffffff,
            ambient: 0x777777,
            specular: 0x333333,
            opacity: 1,
            shininess: shininess,
            shading: shading } );

    DynamicsVisualizer.Materials["lava"] = new THREE.MeshPhongMaterial(
        {   map: lavaTexture,
            bumpMap: lavaTexture,
            bumpScale: 0.5,
            color: 0xffffff,
            ambient: 0x777777,
            specular: 0x333333,
            shininess: shininess,
            shading: shading } );

    DynamicsVisualizer.Materials["moon"] = new THREE.MeshPhongMaterial(
        {   map: moonTexture,
            bumpMap: moonTexture,
            bumpScale: 0.0001,
            color: 0xffffff,
            ambient: 0x777777,
            specular: 0x333333,
            shininess: 0,
            shading: shading } );

    DynamicsVisualizer.Materials["earth"] = new THREE.MeshPhongMaterial(
        {   map: earthTexture,
            color: 0xffffff,
            ambient: 0x777777,
            specular: 0x333333,
            shininess: shininess,
            shading: shading } );

}
loadMaterials();
