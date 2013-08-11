#This module is exclusively for testing Socket connections.

import socket
#from pydy_viz.server import create_socket_server
from server import Server

class TestSocketServer(object):
    def __init__(self):
        self.host = 'localhost'
        self.port = 8000
        self.server = Server(self.port)

    def test_req1(self):
        #We connect to a client socket, and
        #send some GET request to our server socket
        self.connection = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.connection.connect((self.host, self.port))
        self.connection.send('GET /')        
        data = ''
        while data:
            data += self.server.socket.accept().recv(1024)
            
        assert data == open('testindex.html').read()

    def test_req2(self):
        self.connection = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.connection.connect((self.host, self.port))
        self.connection.send('GET /testfile.js')
        data = ''
        while data:
            data += self.server.socket.accept().recv(1024)

        assert  data == open('testfile.js').read()

    def test_req3(self):
        #Test for json file
        pass


