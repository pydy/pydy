#!/usr/bin/env python

from __future__ import print_function

import os
import signal
import socket
import webbrowser
import BaseHTTPServer
from SimpleHTTPServer import SimpleHTTPRequestHandler


class StoppableHTTPServer(BaseHTTPServer.HTTPServer):

    def server_bind(self):
        BaseHTTPServer.HTTPServer.server_bind(self)
        self.socket.settimeout(1)
        self.run = True

    def get_request(self):
        while self.run:
            try:
                sock, addr = self.socket.accept()
                sock.settimeout(None)
                return (sock, addr)
            except socket.timeout:
                pass

    def stop(self):
        self.run = False

    def serve(self):
        while self.run:
            self.handle_request()


class Server(object):
    """
    Parameters
    ----------
    port : integer
        Defines the port on which the server will run.
    scene_file : name of the scene_file generated for visualization
        Used here to display the url
    directory : path of the directory which contains static and scene files.
        Server is started in this directory itself.

    """
    def __init__(self, port=8000, scene_file="Null", directory="static/"):
        self.port = port
        self.scene_file = scene_file
        self.directory = directory

    def run_server(self):
        #change dir to static first.
        os.chdir(self.directory)
        HandlerClass = SimpleHTTPRequestHandler
        ServerClass  = StoppableHTTPServer
        Protocol     = "HTTP/1.0"
        server_address = ('127.0.0.1', self.port)
        HandlerClass.protocol_version = Protocol
        self.httpd = ServerClass(server_address, HandlerClass)
        sa = self.httpd.socket.getsockname()
        print("Serving HTTP on", sa[0], "port", sa[1], "...")
        print("To view visualization, open:\n")
        url = "http://localhost:"+str(sa[1]) + "/index.html?load=" + \
              self.scene_file
        print(url)
        webbrowser.open(url)
        print("hit ctrl+c to stop the server..")
        signal.signal(signal.SIGINT, self.stop_server)
        self.httpd.serve()

    def stop_server(self, signal, frame):
        res = raw_input("Do you really want to exit ([y]/n)? ")
        if res is (None or 'y'):
            self.httpd.stop()
