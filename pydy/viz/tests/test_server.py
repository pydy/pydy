import io
import os
import sys
import signal

# For python 2 and python 3 compatibility
if sys.version_info < (3, 0):
    from SimpleHTTPServer import SimpleHTTPRequestHandler
else:
    from http.server import SimpleHTTPRequestHandler
    raw_input = input

from ..server import StoppableHTTPServer, Server


class TestStoppableHttpServer(object):

    def __init__(self):
        self.request_handler = SimpleHTTPRequestHandler
        self.stoppable_http_server = StoppableHTTPServer(("127.0.0.1", 8080),
                                                         self.request_handler,
                                                         bind_and_activate=False)

    def test_stoppable_http_server(self):
        self.stoppable_http_server.server_bind()
        assert self.stoppable_http_server.running

        self.stoppable_http_server.stop()
        assert not self.stoppable_http_server.running


class TestServer(object):

    def __init__(self):
        directory = os.path.join(os.path.dirname(os.path.abspath(
                                            __file__)), os.pardir, 'static')
        self.test_server = Server(directory=directory,
                            scene_file="js/tests/sample_data/scene_desc.json")

    def test_run_server(self):
        self.test_server.run_server(headless=True)
        assert self.test_server.httpd.running

        sys.stdin = io.StringIO(u'y')
        os.kill(os.getpid(), signal.SIGINT)
        assert not self.test_server.httpd.running
        assert not self.test_server._thread.is_alive()
