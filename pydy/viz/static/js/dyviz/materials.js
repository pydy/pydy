DynamicsVisualizer.Geometries = ['Box', 'Cube', 'Cylinder',
                                 'Cone', 'Sphere', 'Circle',
                                 'Plane', 'Tetrahedron',
                                 'Octahedron', 'Icosahedron',
                                 'Torus', 'TorusKnot']

DynamicsVisualizer.MaterialsList = ['default', 'checkerboard', 'metal',
                                    'dirt', 'foil', 'water', 'grass',
                                    'lava', 'moon', 'earth']

DynamicsVisualizer.Materials = Object.extend(DynamicsVisualizer, {

    getMaterial: function(material) {

        var self = this;

        // TODO : This `pydy-resources` directory should not be hard coded here!
        self.texture_base_path = "pydy-resources/textures/"
        // If IPython,
        if(typeof IPython == "undefined"){
            self.texture_base_path = "textures/";
        }
        self.shininess = 50;
        self.specular = 0x333333;
        self.shading = THREE.SmoothShading;

        switch(material){

            case "default":
                return new THREE.MeshLambertMaterial();

            case "checkerboard":

                self.checkerBoardTexture = THREE.ImageUtils.loadTexture( self.texture_base_path + "checkerboard.jpg" );
                self.checkerBoardTexture.wrapS = self.checkerBoardTexture.wrapT = THREE.RepeatWrapping;
                self.checkerBoardTexture.anisotropy = 4;

                return new THREE.MeshPhongMaterial({
                    map: self.checkerBoardTexture,
                    bumpMap: null,
                    bumpScale: 0,
                    color: 0xffffff,
                    ambient: 0x555555,
                    specular: 0x222222,
                    shininess: 10,
                    shading: self.shading } );

            case "metal":

                self.metalTexture = THREE.ImageUtils.loadTexture( self.texture_base_path + "metal.jpg" );
                self.metalTexture.wrapS = self.metalTexture.wrapT = THREE.RepeatWrapping;
                self.metalTexture.repeat.set(2,2);
                self.metalTexture.anisotropy = 4;

                return new THREE.MeshPhongMaterial({
                    map: self.metalTexture,
                    bumpMap: self.metalTexture,
                    bumpScale: 0.1,
                    color: 0xffffff,
                    ambient: 0x555555,
                    specular: 0x222222,
                    shininess: 70,
                    shading: self.shading } );

            case "dirt":

                self.dirtTexture = THREE.ImageUtils.loadTexture( self.texture_base_path + "dirt.jpg" );
                self.dirtTexture.wrapS = self.dirtTexture.wrapT = THREE.RepeatWrapping;
                self.dirtTexture.repeat.set(1,1);
                self.dirtTexture.anisotropy = 4;

                return new THREE.MeshPhongMaterial({
                    map: self.dirtTexture,
                    bumpMap: self.dirtTexture,
                    bumpScale: 0.05,
                    color: 0xffffff,
                    ambient: 0x555555,
                    specular: 0x222222,
                    shininess: 10,
                    shading: self.shading } );

            case "foil":

                self.waterTexture = THREE.ImageUtils.loadTexture( self.texture_base_path + "water.jpg" );
                self.waterTexture.wrapS = self.waterTexture.wrapT = THREE.RepeatWrapping;
                self.waterTexture.repeat.set(1,1);
                self.waterTexture.anisotropy = 4;

                return new new THREE.MeshPhongMaterial({
                    map: self.waterTexture,
                    bumpMap: self.waterTexture,
                    bumpScale: 0.3,
                    color: 0xccccaa,
                    ambient: 0x555544,
                    specular: 0x777777,
                    shininess: 100,
                    shading: self.shading } );

            case "water":

                self.waterTexture = THREE.ImageUtils.loadTexture( self.texture_base_path + "water.jpg" );
                self.waterTexture.wrapS = self.waterTexture.wrapT = THREE.RepeatWrapping;
                self.waterTexture.repeat.set(1,1);
                self.waterTexture.anisotropy = 4;

                return new THREE.MeshPhongMaterial({
                    map: self.waterTexture,
                    bumpMap: self.waterTexture,
                    bumpScale: 0.05,
                    color: 0x3333aa,
                    ambient: 0x335577,
                    specular: 0x555555,
                    shininess: self.shininess,
                    shading: self.shading } );

            case "grass":

                self.grassTexture = THREE.ImageUtils.loadTexture( self.texture_base_path + "grass.jpg" );
                self.grassTexture.wrapS = self.grassTexture.wrapT = THREE.RepeatWrapping;
                self.grassTexture.repeat.set(1,1);
                self.grassTexture.anisotropy = 4;

                return new THREE.MeshPhongMaterial({
                    map: self.grassTexture,
                    bumpMap: self.grassTexture,
                    bumpScale: 0.05,
                    color: 0xffffff,
                    ambient: 0x777777,
                    specular: 0x333333,
                    opacity: 1,
                    shininess: self.shininess,
                    shading: self.shading } );

            case "lava":

                self.lavaTexture = THREE.ImageUtils.loadTexture( self.texture_base_path + "lavatile.jpg" );
                self.lavaTexture.repeat.set( 4, 2 );
                self.lavaTexture.wrapS = self.lavaTexture.wrapT = THREE.RepeatWrapping;
                self.lavaTexture.anisotropy = 4;

                return new THREE.MeshPhongMaterial({
                    map: self.lavaTexture,
                    bumpMap: self.lavaTexture,
                    bumpScale: 0.5,
                    color: 0xffffff,
                    ambient: 0x777777,
                    specular: 0x333333,
                    shininess: self.shininess,
                    shading: self.shading } );

            case "moon":

                self.moonTexture = THREE.ImageUtils.loadTexture( self.texture_base_path + "moon.jpg" );
                self.moonTexture.wrapS = self.moonTexture.wrapT = THREE.RepeatWrapping;
                self.moonTexture.anisotropy = 4;

                return new THREE.MeshPhongMaterial({
                    map: self.moonTexture,
                    bumpMap: self.moonTexture,
                    bumpScale: 0.0001,
                    color: 0xffffff,
                    ambient: 0x777777,
                    specular: 0x333333,
                    shininess: 0,
                    shading: self.shading } );

            case "earth":

                self.earthTexture = THREE.ImageUtils.loadTexture( self.texture_base_path + "earth.jpg" );
                self.earthTexture.wrapS = self.earthTexture.wrapT = THREE.RepeatWrapping;
                self.earthTexture.anisotropy = 4;

                return new THREE.MeshPhongMaterial({
                    map: self.earthTexture,
                    color: 0xffffff,
                    ambient: 0x777777,
                    specular: 0x333333,
                    shininess: self.shininess,
                    shading: self.shading } );
        }
    }
});
