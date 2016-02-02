from SimpleHTTPServer import SimpleHTTPRequestHandler
from BaseHTTPServer import HTTPServer

from ..server import StoppableHTTPServer


class TestStoppableHttpServer(object):

    def __init__(self):
        self.request_handler = SimpleHTTPRequestHandler
        self.http_server = HTTPServer

        self.stoppable_http_server = StoppableHTTPServer(("127.0.0.1", 8080),
                                                         self.request_handler,
                                                         bind_and_activate=False)

    def test_stoppable_http_server(self):
        self.stoppable_http_server.server_bind()
        assert self.stoppable_http_server.run is True

        self.stoppable_http_server.stop()
        assert self.stoppable_http_server.run is False
