#!/usr/bin/env python

import os
import webbrowser
import BaseHTTPServer
from SimpleHTTPServer import SimpleHTTPRequestHandler


def run_server(port=8000,scene_file="Null"):
    #change dir to static first.
    os.chdir("static/")
    HandlerClass = SimpleHTTPRequestHandler
    ServerClass  = BaseHTTPServer.HTTPServer
    Protocol     = "HTTP/1.0"
    server_address = ('127.0.0.1', port)
    HandlerClass.protocol_version = Protocol
    httpd = ServerClass(server_address, HandlerClass)
    sa = httpd.socket.getsockname()
    print "Serving HTTP on", sa[0], "port", sa[1], "..."
    print "hit ctrl+c to stop the server.."
    print "To view visualization, open:\n"
    url = "http://localhost:"+str(sa[1]) + "/index.html?load="+scene_file
    print url
    webbrowser.open(url)
    httpd.serve_forever()


if __name__ == "__main__":
    run_server()