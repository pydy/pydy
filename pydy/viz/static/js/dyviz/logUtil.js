
(function($) {

DynamicsVisualizer.Logger = Object.extend(DynamicsVisualizer, {
    
    Success: function(content){
        toAppend = '<li class="log-success"' + content + '</li>'
        $("#log-container").append(toAppend);
 
    },

    Info: function(content){

    },
    Warning: function(content){

    },

    Danger: function(content){


    }

})

})(jQuery);

   


