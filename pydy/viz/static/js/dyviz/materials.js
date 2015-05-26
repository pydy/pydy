DynamicsVisualizer.Geometries = ['Cube', 'Cylinder',
                                 'Cone', 'Sphere', 'Circle',
                                 'Plane', 'Tetrahedron',
                                 'Octahedron', 'Icosahedron',
                                 'Torus', 'TorusKnot']

DynamicsVisualizer.Materials = Object.extend(DynamicsVisualizer, {

    loadTexture: function(){

        var self = this;

        self.texture_base_path = "static/textures/"
        // If IPython,
        if(typeof IPython == "undefined"){
            self.texture_base_path = "textures/";
        }

        self.checkerBoardTexture = THREE.ImageUtils.loadTexture( self.texture_base_path + "checkerboard.jpg" );
        self.checkerBoardTexture.wrapS = self.checkerBoardTexture.wrapT = THREE.RepeatWrapping;
        self.checkerBoardTexture.anisotropy = 4;

        self.metalTexture = THREE.ImageUtils.loadTexture( self.texture_base_path + "metal.jpg" );
        self.metalTexture.wrapS = self.metalTexture.wrapT = THREE.RepeatWrapping;
        self.metalTexture.repeat.set(2,2);
        self.metalTexture.anisotropy = 4;

        self.moonTexture = THREE.ImageUtils.loadTexture( self.texture_base_path + "moon.jpg" );
        self.moonTexture.wrapS = self.moonTexture.wrapT = THREE.RepeatWrapping;
        self.moonTexture.anisotropy = 4;

        self.earthTexture = THREE.ImageUtils.loadTexture( self.texture_base_path + "earth.jpg" );
        self.earthTexture.wrapS = self.earthTexture.wrapT = THREE.RepeatWrapping;
        self.earthTexture.anisotropy = 4;

        self.grassTexture = THREE.ImageUtils.loadTexture( self.texture_base_path + "grass.jpg" );
        self.grassTexture.wrapS = self.grassTexture.wrapT = THREE.RepeatWrapping;
        self.grassTexture.repeat.set(1,1);
        self.grassTexture.anisotropy = 4;

        self.dirtTexture = THREE.ImageUtils.loadTexture( self.texture_base_path + "dirt.jpg" );
        self.dirtTexture.wrapS = self.dirtTexture.wrapT = THREE.RepeatWrapping;
        self.dirtTexture.repeat.set(1,1);
        self.dirtTexture.anisotropy = 4;

        self.waterTexture = THREE.ImageUtils.loadTexture( self.texture_base_path + "water.jpg" );
        self.waterTexture.wrapS = self.waterTexture.wrapT = THREE.RepeatWrapping;
        self.waterTexture.repeat.set(1,1);
        self.waterTexture.anisotropy = 4;

        self.lavaTexture = THREE.ImageUtils.loadTexture( self.texture_base_path + "lavatile.jpg" );
        self.lavaTexture.repeat.set( 4, 2 );
        self.lavaTexture.wrapS = self.lavaTexture.wrapT = THREE.RepeatWrapping;
        self.lavaTexture.anisotropy = 4;

        self.shininess = 50;
        self.specular = 0x333333;
        self.shading = THREE.SmoothShading;
    },

    getMaterial: function(material) {

        switch(material){

            case "default":
                return new THREE.MeshLambertMaterial();

            case "checkerboard":
                return new THREE.MeshPhongMaterial({
                    map: checkerBoardTexture,
                    bumpMap: null,
                    bumpScale: 0,
                    color: 0xffffff,
                    ambient: 0x555555,
                    specular: 0x222222,
                    shininess: 10,
                    shading: shading } );

            case "metal":
                return new THREE.MeshPhongMaterial({
                    map: self.metalTexture,
                    bumpMap: self.metalTexture,
                    bumpScale: 0.1,
                    color: 0xffffff,
                    ambient: 0x555555,
                    specular: 0x222222,
                    shininess: 70,
                    shading: shading } );

            case "dirt":
                return new THREE.MeshPhongMaterial({
                    map: self.dirtTexture,
                    bumpMap: self.dirtTexture,
                    bumpScale: 0.05,
                    color: 0xffffff,
                    ambient: 0x555555,
                    specular: 0x222222,
                    shininess: 10,
                    shading: shading } );

            case "foil":
                return new new THREE.MeshPhongMaterial({
                    map: self.waterTexture,
                    bumpMap: self.waterTexture,
                    bumpScale: 0.3,
                    color: 0xccccaa,
                    ambient: 0x555544,
                    specular: 0x777777,
                    shininess: 100,
                    shading: shading } );

            case "water":
                return new THREE.MeshPhongMaterial({
                    map: self.waterTexture,
                    bumpMap: self.waterTexture,
                    bumpScale: 0.05,
                    color: 0x3333aa,
                    ambient: 0x335577,
                    specular: 0x555555,
                    shininess: shininess,
                    shading: shading } );

            case "grass":
                return new THREE.MeshPhongMaterial({
                    map: self.grassTexture,
                    bumpMap: self.grassTexture,
                    bumpScale: 0.05,
                    color: 0xffffff,
                    ambient: 0x777777,
                    specular: 0x333333,
                    opacity: 1,
                    shininess: shininess,
                    shading: shading } );

            case "lava":
                return new THREE.MeshPhongMaterial({
                    map: self.lavaTexture,
                    bumpMap: self.lavaTexture,
                    bumpScale: 0.5,
                    color: 0xffffff,
                    ambient: 0x777777,
                    specular: 0x333333,
                    shininess: shininess,
                    shading: shading } );

            case "moon":
                return new THREE.MeshPhongMaterial({
                    map: self.moonTexture,
                    bumpMap: self.moonTexture,
                    bumpScale: 0.0001,
                    color: 0xffffff,
                    ambient: 0x777777,
                    specular: 0x333333,
                    shininess: 0,
                    shading: shading } );

            case "earth":
                return new THREE.MeshPhongMaterial({
                    map: self.earthTexture,
                    color: 0xffffff,
                    ambient: 0x777777,
                    specular: 0x333333,
                    shininess: shininess,
                    shading: shading } );
        }
    }
});

// Load all the textures for materials
DynamicsVisualizer.Materials.loadTexture();
