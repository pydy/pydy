#This module is exclusively for testing Socket connections.

import socket
from pydy_viz.server import create_socket_server

class TestSocketServer(object):
    def __init__(self):
        self.host = 'localhost'
        self.port = 8000
        create_socket_server(8000)

    def test_req1(self):
        
        #We connect to a client socket, and 
        #send some GET request to our server socket
        
        assert data == #Some index.html data

    def test_req2(self):
        self.connection = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.connection.connect((self.host, self.port))
        connection.send('GET /testfile.js')
        data = s.recv(1024)
        s.close()

        assert  data == open('testfile.js'.read()
    
    def test_req3(self):
        #Test for json file
        pass    


