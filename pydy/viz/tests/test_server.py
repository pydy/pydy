import sys

# For python 2 and python 3 compatibility
if sys.version_info < (3, 0):
    from SimpleHTTPServer import SimpleHTTPRequestHandler
else:
    from http.server import SimpleHTTPRequestHandler
    raw_input = input

from ..server import StoppableHTTPServer


class TestStoppableHttpServer(object):

    def __init__(self):
        self.request_handler = SimpleHTTPRequestHandler
        self.stoppable_http_server = StoppableHTTPServer(("127.0.0.1", 8080),
                                                         self.request_handler,
                                                         bind_and_activate=False)

    def test_stoppable_http_server(self):
        self.stoppable_http_server.server_bind()
        assert self.stoppable_http_server.run

        self.stoppable_http_server.stop()
        assert not self.stoppable_http_server.run
