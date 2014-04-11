import os
import socket
import sys
import threading

__all__ = ['Server']

class Server(threading.Thread):
    """
    A basic Socket server.
    This server is used for fetching
    static files from the pydy_viz
    source and rendering them to
    the browser.

    """
    def __init__(self, json='data.json'):
        """
        Initiate a server instance.

        Parameters
        ==========
        json : str, optional
        path to the saved json file
        for visualization, relative to
        current working directory.


        """
        threading.Thread.__init__(self)
        self.saved_json_file = json
        self.port = 8000
        host = ''
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            self.socket.bind((host, self.port))
        except:
            self.port += 1
            try:
                self.socket.bind((host, self.port))
            except:
                self.port += 1
                self.socket.bind((host, self.port))
        self.socket.listen(5)

    def _parse_data(self, data):

        static_path = os.path.dirname(__file__)
        static_path = os.path.join(static_path, 'static')
        try:
            request = data.split(' ')[1]
        except IndexError:
			#If error occurs in parsing a request,
			#Better to reload the page, to avoid broken
			#javascripts
			request = '/'
        if request == '/':
        #If requested for http://localhost:port/
        #Send index.html file
            file_path = os.path.join(static_path, 'index.html')
            send_buffer = ''

        elif request == '/data.json':
        #If data.json is requested, get it from scene method
            file_path = os.path.join(os.getcwd(), self.saved_json_file)
            send_buffer = 'var JSONObj = '

        elif request == '/close-server':
            print "Server closed successfully!"
            self.close()


        else:
        #Else, try to use relative path from url for other files
            file_path_list = request.split('/')
            for val in file_path_list:
                static_path = os.path.join(static_path, val)
            file_path = static_path
            send_buffer = ''

        try:
            send_buffer += open(file_path).read()

        except IOError:
            pass
            #print '''404 File not found. Sent No Data'''

        return send_buffer

    def listen_once(self):
        conn, addr = self.socket.accept()
        data = conn.recv(1024)
        sent_data = self._parse_data(data)
        conn.send(sent_data)
        return sent_data

    def run(self):
        print 'server started successfully, on port:', self.port

        while 1:
            try:
                self.socket.listen(1)
                conn, addr = self.socket.accept()
                data = conn.recv(1024)
                sent_data = self._parse_data(data)
                conn.send('HTTP/1.1 200 OK\r\n\r\n')
                conn.send(sent_data)
                conn.close()
            except KeyboardInterrupt:
                print "Are you sure you want to shutdown[Y/N]?"
                a = raw_input()
                if a == "Y" or a == "y":
                    self.close()
                else:
                    pass
            except socket.error, e:
                if isinstance(e.args, tuple):
                    print "errno is %d" % e[0]
                if e[0] == errno.EPIPE:
                    # remote peer disconnected
                    print "Detected remote disconnect"

    def close(self):
        self.socket.shutdown(socket.SHUT_RDWR)
        sys.exit()

if __name__ == "__main__":
    a = Server()
    a.run()
