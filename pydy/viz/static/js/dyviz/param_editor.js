DynamicsVisualizer.ParamEditor = Object.extend(DynamicsVisualizer, {
	
	openDialog: function(id){

		/** 
		  * 
		**/ 
        var self = this;

        var toLoad = self._scene.getObjectByName(parseInt(id));
        toLoad = toLoad["object-info"];
        // Remove the old blinker and reset material to default..
        window.clearInterval(self.blinkId);
        //self._blinker.material = self._old_material;
        //
        self.Scene._blink(parseInt(id));
        var mainDiv = jQuery('<div/>',{id: "object-"+ toLoad.simulation_id, style: 'display:none;'}); 
    	
    	// for name..
        var div1 = jQuery('<div />',{class: 'input-group'});
    	div1.append('<span class="input-group-addon">Name</span>');
    	div1.append(jQuery('<input />',{ type:'text', id: "_name", class: 'form-control', value: toLoad.name}));

        // for color..
    	div1.append('<span class="input-group-addon">Color</span>');
    	div1.append(jQuery('<input />',{ type:'text', id: "_color", class: 'form-control', value: toLoad.color}));
        
        // for material.. a dropdown..

        var div_material = jQuery('<select />',{class: 'form-control', id:"_material"});
        div_material.append('<option value="' + toLoad.material + '">' + toLoad.material + '</option>');
        for(var i in self.Materials){
        	div_material.append('<option value="' + i +  '">' + i + '</option>');
        }

        var div_geom = jQuery('<select />',{class: 'form-control', id:"_geometry"});
        div_geom.append('<option value="' + toLoad.type + '">' + toLoad.type + '</option>');
        for(var i=0;i<self.Geometries.length; i++){
        	div_geom.append('<option value="' + self.Geometries[i] +  '">' + self.Geometries[i] + '</option>');
        }

        console.log(jQuery("#geometry"));
        
        //rest geom. params depend on shape..
        var div2 = jQuery('<div />',{class: 'input-group', id: "geom-params"});
        
        mainDiv.append(div1);
        mainDiv.append('<hr/><span class="input-group-addon">Material</span>');
        mainDiv.append(div_material);
        mainDiv.append('<hr/><span class="input-group-addon">Geometry</span>');
        mainDiv.append(div_geom);
        mainDiv.append("<hr />");
        mainDiv.append(div2);
        // finally a button..
        mainDiv.append('<hr /><button id="apply-' + id +  '" class="btn btn-primary btn-small">Apply</button>');
        jQuery("#object-dialog").html(mainDiv);

        // show after whole div is populated..
        mainDiv.fadeIn("slow");

        jQuery("#_geometry").change( function(){
        
        	self.ParamEditor._addGeometryFor(jQuery(this).val());
        });

        // finally activate button..
        jQuery("#apply-" + id).click(function(){
			self.ParamEditor.applySceneInfo(jQuery(this).attr("id").split("-").slice(-1)[0]);
			
		});

        self.ParamEditor._addGeometryFor(toLoad);
        jQuery("#close-object-dialog").removeClass("disabled");


    },

    applySceneInfo: function(id){
        /** 
          * To apply scene info, we modify the object and then 
          * call Scene._addIndividualObject on it
          * 
        **/
    	var self = this;
        // remove blinker..
        window.clearInterval(self.blinkId);
        var int_id = parseInt(id);
    	var updated_object = {};

        updated_object.name = jQuery("#_name").val();
        updated_object.color = jQuery("#_color").val();
        updated_object.material = jQuery("#_material").val();
        updated_object.type = jQuery("#_geometry").val();
        updated_object.simulation_id = int_id;
        

        //Now by type, use switch case..
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

            // TODO for Mesh..     


        }

        console.log("updated object: ");
        console.log(updated_object);
        
        jQuery.extend(true,self.model.objects[int_id],updated_object);
        console.log(self.model)
        self.Scene._removeAll();
        self.Scene.addObjects();
        self.Scene.addCameras();
        self.Scene.addLights();
        self.loadUIElements();
        //self.Scene._addIndividualObject(updated_object);


    },

    _addGeometryFor: function(toLoad){
    	var self = this;
        var div2 = jQuery("#geom-params");
        div2.html(" ");
        
        switch(toLoad.type || toLoad){
        	// TODO all objects..
        	case "Sphere":
            case "Circle":
            case "Tetrahedron":
            case "Octahedron":
            case "Icosahedron":
        	    div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Radius'));
    	        div2.append(jQuery('<input />',{ type:'text', id: "_radius", class: 'form-control', value: toLoad.radius || 0.0}));
        	    break;

        	case "Cylinder":
            case "Cone":
        	    div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Radius'));
    	        div2.append(jQuery('<input />',{ type:'text', id: "_radius", class: 'form-control', value: toLoad.radius  || 0.0}));

    	        div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Length'));
    	        div2.append(jQuery('<input />',{ type:'text', id: "_length", class: 'form-control', value: toLoad.length || 0.0}));
        	    break;

            case "Torus":
            case "TorusKnot":
                div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Radius'));
                div2.append(jQuery('<input />',{ type:'text', id: "_radius", class: 'form-control', value: toLoad.radius  || 0.0}));

                div2.append(jQuery('<span \>',{ class:'input-group-addon',}).html('Tube Radius'));
                div2.append(jQuery('<input />',{ type:'text', id: "_tubeRadius", class: 'form-control', value: toLoad.tube_radius || 0.0}));
                break;    


        }
    	
    },

    showModel: function(){
        // update editor..

        self.editor.getDoc().setValue(JSON.stringify(self.model,null,4));
        jQuery("#model-loader-wrapper").slideIn();

        self.editor.refresh();
    }

});