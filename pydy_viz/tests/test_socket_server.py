#This module is exclusively for testing Socket connections.

import socket
from pydy_viz.server import Server
import pydy_viz
import os

class TestSocketServer(object):
    def __init__(self):
        self.host = 'localhost'
        self.server = Server()

    def test_req1(self):
        #We connect to a client socket, and
        #send some GET request to our server socket
        self.connection = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.connection.connect((self.host, self.server.port))
        self.connection.send('GET /')

        data = self.server.listen_once()
        path = os.path.join(os.path.dirname(pydy_viz.__file__), 'static', 'index.html')
        print data
        assert data == open(path).read()

    def test_req2(self):
        self.connection = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.connection.connect((self.host, self.server.port))
        self.connection.send('GET /js/canvas/initialize.js')

        data = self.server.listen_once()
        path = os.path.join(os.path.dirname(pydy_viz.__file__), 'static', 'js', 'canvas', 'initialize.js')
        print data
        assert data == open(path).read()


