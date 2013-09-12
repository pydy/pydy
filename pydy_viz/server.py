__all__ = ['Server']

import socket, threading, os
import pydy_viz
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
            self.port+=1
            try:
                self.socket.bind((host, self.port))
            except:
                self.port+=1
                self.socket.bind((host, self.port))
        self.socket.listen(5)

    def _parse_data(self, data):

        static_path = os.path.dirname(pydy_viz.__file__)
        static_path = os.path.join(static_path, 'static')
        print "static path", static_path
        request = data.split(' ')[1]
        print 'request:%s'% request
        if request == '/':
        #If requested to http://localhost:port/
        #Send index.html file
            file_path = os.path.join(static_path, 'index.html')
            send_buffer = ''

        elif request == '/data.json':
            print 'data file requested'
        #If data.json is requested, get it from scene method
            file_path = os.path.join(os.getcwd(), self.saved_json_file)
            send_buffer = 'var JSONObj = '

        else:
        #Else, try to use relative path from url for other files
            file_path_list = request.split('/')
            for val in file_path_list:
                static_path = os.path.join(static_path, val)
            file_path = static_path
            print "file path : " + file_path
            send_buffer = ''

        try:
            send_buffer += open(file_path).read()

        except IOError:
            print '''404 File not found. Sent No Data'''


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
                self.socket.listen(1)
                conn, addr = self.socket.accept()
                self.data = conn.recv(1024)
                sent_data = self._parse_data(self.data)
                conn.send(sent_data)
                print 'sent data'
                conn.close()

if __name__ == "__main__":
    a = Server()
    a.run()
