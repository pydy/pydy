import sys
import subprocess
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
        assert self.stoppable_http_server.run

        self.stoppable_http_server.stop()
        assert not self.stoppable_http_server.run


class TestServer(object):

    def __init__(self):
        self.test_server = Server(directory="../static/",
                scene_file="js/tests/sample_data/scene_desc.json")

    def test_run_server(self):
        assert self.test_server.directory == "../static/"
        assert self.test_server.port == 8000

        proc = subprocess.Popen(self.test_server.run_server(),
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)

        (stdout, stderr) = proc.communicate()
        sa = self.test_server.httpd.socket.getsockname()
        assert sa[0] == "127.0.0.1"
        assert self.test_server.httpd.run

        proc.communicate(signal.SIGINT, timeout=1)
        (stdout, stderr) = proc.communicate()

        proc.communicate("y", timeout=1)
        (stdout, stderr) = proc.communicate()

        proc.terminate()

        assert not self.test_server.httpd.run
