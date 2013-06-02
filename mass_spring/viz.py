from IPython.display import Javascript, display
import re
import os


# Load required assets on import
path = os.path.dirname(os.path.realpath(__file__))
js = {"three": "js/three.min.js",
      "controls": "js/TrackballControls.js",
      "bob": "js/bob.js"}
for key, filename in js.items():
    with open(os.path.join(path, filename)) as in_js:
        js[key] = in_js.read()

# This is the only way I found to use local copies of js libraries in IPython
js["script"] = js["three"] + js["controls"] + js["bob"]


def visualize(coordinates):
	for coordinate in coordinates:
		
		
		
		replacer = {'y':str(coordinate[0]),'time':str(coordinate[1])}
		
		js_data = re.sub("#\(\w+\)", lambda m: replacer[m.group()[2:-1]], js["script"])
		
		display(Javascript(js_data))
	
