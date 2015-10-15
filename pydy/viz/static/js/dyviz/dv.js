// Dynamics Visualizer main class.

var DynamicsVisualizer = {};

DynamicsVisualizer = Class.create({
    /**
      * DV is the main class for Dynamics Visualizer.
      * it contains methods to set up a default UI, and
      * maps buttons' `onClick` to functions.
    **/

    _initialize: function(){
        /**
          * Checks whether the browser supports webGLs, and
          * initializes the DynamicVisualizer object.
        **/
        var self = this;
        console.log("[PyDy INFO]: initializing Visualizer");
        if(!self.isWebGLCompatible()){
            console.log("[PyDy ALERT]: Incompatible browser!");
            alert("The browser does not seems to be webgl compatible! " +
                "Please check here for browser compatibility: http://caniuse.com/webgl ");
            return false;
        }

        var sceneFileURI = self.getQueryString("load") || jQuery("#json-input").val();
        if(sceneFileURI){
            console.log("[PyDy INFO]: Found scene desc from URL");
            jQuery("#json-input").val(sceneFileURI);
            self.sceneFilePath = sceneFileURI;
            console.log("[PyDy INFO]: Loading scene JSON file:" + self.sceneFilePath);
            self.Parser.loadScene();
        }

    },


    isWebGLCompatible: function(){
        /**
          * Checks whether the browser used is
          * compatible for handling webGL based
          * animations.
          * Requires external script: Modernizr.js
          *
        **/

        if (!Modernizr.canvas || !Modernizr.webgl) return false;
        else return true;
    },

    activateUIControls: function(){
        /**
          * This method adds functions to the UI buttons
          * It should be **strictly** called after the
          * other DynamicsVisualizer sub-modules are loaded
          * in the browser, else certain functionality will
          * be(not might be!) hindered.
        **/


        var self = this;
        jQuery("#simulation-load").click(function(){
            self.sceneFilePath = jQuery("#json-input").val();
            console.log("[PyDy INFO]: Loading scene JSON file:" + self.sceneFilePath);
            self.Parser.loadScene();
        });


        self._slider = jQuery("#time-slider").slider({min:0, max:100, step:1, handle:"square", value:0});
        self._slider.on('slide',function(ev) {
            var val = ev.value;
            var len = self._timeArray.length;
            var i = 0;
            var gotValue = false;
            for( i=0;i<self._timeArray.length && !gotValue;i++){
                var percent = (self._timeArray[i]/self._timeArray[len-1])*100;
                if(val <= percent){ gotValue = true; break; }
            }
            self.currentTime = self._timeArray[i];
            self.Scene.setAnimationTime(self._timeArray[i]);
        });

        jQuery("#resetControls").click(function(){
            self.scene._resetControls();
        });

        jQuery("#play-animation").click(function(){
            self.Scene.runAnimation();

        });
        jQuery("#pause-animation").click(function(){
            self.Scene.pauseAnimation();

        });
        jQuery("#stop-animation").click(function(){
            self.Scene.stopAnimation();

        });
        jQuery("#close-object-dialog").click(function(){
            jQuery("#object-dialog").html(" ");
            jQuery(this).addClass("disabled");
            jQuery("#object-dialog").css("display", "none");
        });

        jQuery("#show-model").click(function(){
            jQuery("#model-loader-wrapper").slideDown();
            self.editor.refresh();
            // Make JSON downloadable..
            var data = "text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(self.model,null,"    "));
            jQuery('<a href="data:' + data +
                '" download="scene_desc.json" class="btn btn-success btn-large"> \
                <i class="icon-white icon-download-alt">download JSON</a>').appendTo('#download-json');
            jQuery(this).addClass("disabled");
        });

        jQuery("#close-model-dialog").click(function(){
            jQuery("#model-loader-wrapper").slideUp();
            jQuery("#download-json").html(" ");
            jQuery("#show-model").removeClass("disabled")

        });
        console.log("[PyDy INFO]: Activated UI controls");

    },

    loadUIElements: function(){
        /**
          * This method loads UI elements
          * which can be loaded only **after**
          * scene JSON is loaded onto canvas.
        **/
        var self = this;
        jQuery("#play-animation").removeClass("disabled");
        jQuery("#show-model").removeClass("disabled");
        var objs = self.model.objects;

        jQuery("#object-dropdown").find("li").remove(); // clean old dropdown list.

        for(var obj in objs){
            var toAppend = '<li><a id="'+ objs[obj].simulation_id +
                           '" href="#">' + objs[obj].name + '</a></li>';
            jQuery("#object-dropdown").append(toAppend);
            // adding click functions to all dropdown objs.
            jQuery("#" + objs[obj].simulation_id).click(function(){
                self.ParamEditor.openDialog(jQuery(this).attr("id"));
                jQuery("#object-dialog").css("display", "block");
            });
        }

        var constants = self.model.constant_map;
        var div = jQuery("#simulation-params").fadeOut();
        div.html(" "); // clear html first

        for(var i in constants){
            div.append('<span class="input-group-addon">' + i + '</span>');
            div.append(jQuery('<input />',{ type:'text', id: i, class: 'form-control', value: constants[i]}));
        }
        div.fadeIn();

        // enable CodeMirror...
        jQuery("#model-loader").html(" ");
        if(!self.editor){
            self.editor = CodeMirror.fromTextArea(document.getElementById('model-loader'), {
                height: "30em",
                mode: {name: "javascript", json: true},
                theme: "base16-light",
                textWrapping: true
            });
        }
        self.editor.getDoc().setValue(JSON.stringify(self.model,null,4));

        // Get animation Speed..
        self.animSpeed = jQuery("#anim-speed").val();
    },


    getBasePath: function(){
        /**
          * Returns the base path of
          * the loaded Scene file.
        **/
        var self = this;

        var slashes_fixed = self.sceneFilePath.replace(/\\/g, "/");
        return slashes_fixed.split("/").slice(0,-1).join("/") + "/";
    },

    getFileExtenstion: function(){
        /**
          * Returns the extension of
          * the uploaded Scene file.
        **/
        var self = this;
        return self.sceneFilePath.split(".").slice(-1)[0].toLowerCase();

    },

   getQueryString: function(key){
        /**
          * Returns the GET Parameter from url corresponding
          * to `key`
        **/
        key = key.replace(/[\[]/, "\\[").replace(/[\]]/, "\\]");
        var regex = new RegExp("[\\?&]" + key + "=([^&#]*)"),
        results = regex.exec(location.search);
        return results == null ? "" : decodeURIComponent(results[1].replace(/\+/g, " "));
    }

});

var DynamicsVisualizer = new DynamicsVisualizer();
