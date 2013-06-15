from IPython.display import Javascript, display, HTML
import re
import os


# Get to the current directory
path = os.path.dirname(os.path.realpath(__file__))

#Loading basic Js libraries
display(HTML("""<script src="files/js/three.min.js"></script></script> 
                <script src="files/js/TrackballControls.js"></script>
                """))

display(Javascript('console.log("Imported JS libraries");'))

#opening the Js file
filename  = 'js/bob.js'      
js  = open(os.path.join(path, filename)).read()
      



def visualize(coordinate_time):
    #we will replace the data in the js_data file
    js_data = re.sub("#\(coordinates_time\)", str(coordinate_time),  js)
    display(Javascript(js_data))
