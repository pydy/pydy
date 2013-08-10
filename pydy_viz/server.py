import sys
import BaseHTTPServer
from SimpleHTTPServer import SimpleHTTPRequestHandler


def create_server(port=8000):
    HandlerClass = SimpleHTTPRequestHandler
    ServerClass  = BaseHTTPServer.HTTPServer
    Protocol     = "HTTP/1.0"

    server_address = ('127.0.0.1', port)

    HandlerClass.protocol_version = Protocol
    httpd = ServerClass(server_address, HandlerClass)

    sa = httpd.socket.getsockname()
    print 'Your server is running on http://127.0.0.1:8000/ \n \
           Visit in browser to see your visulization in all its glory'
    httpd.serve_forever()
