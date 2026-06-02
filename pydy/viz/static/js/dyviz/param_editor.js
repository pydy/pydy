DynamicsVisualizer.ParamEditor = Object.extend(DynamicsVisualizer, {

    openDialog: function(id){

        /**
          * This function takes object's id
          * as the argument, and populates the
          * edit objects dialog box.
        **/
        var self = this;

        var toLoad = self._scene.getObjectByName(parseInt(id));
        toLoad = toLoad["object-info"];
        window.clearInterval(self.blinkId);
        if(typeof self._blinker != "undefined"){
            self._blinker.visible = true;
        }
        self.Scene._blink(parseInt(id));
        var mainDiv = jQuery('<div/>',{id: "object-"+ toLoad.simulation_id, style: 'display:none;'});

        var div1 = jQuery('<div />',{class: 'input-group'});
        div1.append('<span class="input-group-addon">Name</span>');
        div1.append(jQuery('<input />',{ type:'text', id: "_name", class: 'form-control', value: toLoad.name}));

        div1.append('<span class="input-group-addon">Color</span>');
        div1.append(jQuery('<input />',{ type:'text', id: "_color", class: 'form-control', value: toLoad.color}));

        var div_material = jQuery('<select />',{class: 'form-control', id:"_material"});
        for(var i=0;i<self.MaterialsList.length; i++){
            if(self.MaterialsList[i] == toLoad.type){
                div_material.append('<option value="' + self.MaterialsList[i] +  '" selected="selected">' + self.MaterialsList[i] + '</option>');
            }  else {
                div_material.append('<option value="' + self.MaterialsList[i] +  '">' + self.MaterialsList[i] + '</option>');
            }
        }

        var div_geom = jQuery('<select />',{class: 'form-control', id:"_geometry"});
        for(var i=0;i<self.Geometries.length; i++){
            if(self.Geometries[i] == toLoad.type){
                div_geom.append('<option value="' + self.Geometries[i] +  '" selected="selected">' + self.Geometries[i] + '</option>');
            }  else {
                div_geom.append('<option value="' + self.Geometries[i] +  '">' + self.Geometries[i] + '</option>');
            }
        }

        var div2 = jQuery('<div />',{class: 'input-group', id: "geom-params"});

        mainDiv.append(div1);
        mainDiv.append('<hr/><span class="input-group-addon">Material</span>');
        mainDiv.append(div_material);
        mainDiv.append('<hr/><span class="input-group-addon">Geometry</span>');
        mainDiv.append(div_geom);
        mainDiv.append("<hr />");
        mainDiv.append(div2);

        mainDiv.append('<hr /><button id="apply-' + id +  '" class="btn btn-primary btn-small">Apply</button>');
        jQuery("#object-dialog").html(mainDiv);

        mainDiv.fadeIn("slow");
        jQuery("#_geometry").change( function(){
            self.ParamEditor._addGeometryFor(jQuery(this).val());
        });

        jQuery("#apply-" + id).click(function(){
            self.ParamEditor.applySceneInfo(jQuery(this).attr("id").split("-").slice(-1)[0]);
        });

        self.ParamEditor._addGeometryFor(toLoad);
        jQuery("#close-object-dialog").removeClass("disabled");
    },

    applySceneInfo: function(id){
        /**
          * This object applies the changes made in
          * the edit objects dialog box to self.model
          * and then renders the model onto canvas.
          * It takes the id of the object as its argument.
        **/
        var self = this;
        window.clearInterval(self.blinkId);
        self._blinker.visible = true;

        var int_id = parseInt(id);
        var updated_object = {};

        updated_object.name = jQuery("#_name").val();
        updated_object.color = jQuery("#_color").val();
        updated_object.material = jQuery("#_material").val();
        updated_object.type = jQuery("#_geometry").val();
        updated_object.simulation_id = int_id;

        switch(updated_object.type){
            case "Sphere":
            case "Circle":
            case "Tetrahedron":
            case "Octahedron":
            case "Icosahedron":
                updated_object.radius = jQuery("#_radius").val()
                break;

            case "Cylinder":
            case "Cone":
                updated_object.radius = jQuery("#_radius").val()
                updated_object.length = jQuery("#_length").val()
                break;

            case "Torus":
            case "TorusKnot":
                updated_object.radius = jQuery("#_radius").val()
                updated_object.tube_radius = jQuery("#_tubeRadius").val()
                break;

            case "Plane":
                updated_object.width = jQuery("#_width").val()
                updated_object.length = jQuery("#_length").val()
                break;

            case "Cube":
                updated_object.length = jQuery("#_length").val()
                break;

            case "Box":
                updated_object.width = jQuery("#_width").val()
                updated_object.height = jQuery("#_height").val()
                updated_object.depth = jQuery("#_depth").val()
                break;

        }

        jQuery.extend(true,self.model.objects[int_id],updated_object);
        self.Scene.addObjects();
        self.Scene.addLights();
        self.loadUIElements();
    },

    _addGeometryFor: function(toLoad){
        /**
          * Adds geometry info for a particular
          * object onto the edit objects dialog
          * box. Takes the object as the argument.
        **/
        var self = this;
        var div2 = jQuery("#geom-params");
        div2.html(" ");

        switch(toLoad.type || toLoad){
            case "Sphere":
            case "Circle":
            case "Tetrahedron":
            case "Octahedron":
            case "Icosahedron":
                div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Radius'));
                div2.append(jQuery('<input />',{ type:'text', id: "_radius", class: 'form-control', value: toLoad.radius || 1.0}));
                break;

            case "Cylinder":
            case "Cone":
                div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Radius'));
                div2.append(jQuery('<input />',{ type:'text', id: "_radius", class: 'form-control', value: toLoad.radius  || 1.0}));

                div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Length'));
                div2.append(jQuery('<input />',{ type:'text', id: "_length", class: 'form-control', value: toLoad.length || 1.0}));
                break;

            case "Torus":
            case "TorusKnot":
                div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Radius'));
                div2.append(jQuery('<input />',{ type:'text', id: "_radius", class: 'form-control', value: toLoad.radius  || 1.0}));

                div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Tube Radius'));
                div2.append(jQuery('<input />',{ type:'text', id: "_tubeRadius", class: 'form-control', value: toLoad.tube_radius || 1.0}));
                break;

            case "Plane":
                div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Width'));
                div2.append(jQuery('<input />',{ type:'text', id: "_width", class: 'form-control', value: toLoad.width  || 1.0}));

                div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Length'));
                div2.append(jQuery('<input />',{ type:'text', id: "_length", class: 'form-control', value: toLoad.length || 1.0}));
                break;

            case "Cube":
                div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Length'));
                div2.append(jQuery('<input />',{ type:'text', id: "_length", class: 'form-control', value: toLoad.length || 1.0}));
                break;

            case "Box":
                div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Width'));
                div2.append(jQuery('<input />',{ type:'text', id: "_width", class: 'form-control', value: toLoad.width || 1.0}));
                
                div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Height'));
                div2.append(jQuery('<input />',{ type:'text', id: "_height", class: 'form-control', value: toLoad.height || 0.5}));
                
                div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Depth'));
                div2.append(jQuery('<input />',{ type:'text', id: "_depth", class: 'form-control', value: toLoad.depth || 0.5}));
                break;


        }
    },

    showModel: function(){
        /**
          * Updates the codemirror instance with
          * the updated model, and shows it in the
          * UI.
         **/
        self.editor.getDoc().setValue(JSON.stringify(self.model,null,4));
        jQuery("#model-loader-wrapper").slideIn();

        self.editor.refresh();
    }
});
