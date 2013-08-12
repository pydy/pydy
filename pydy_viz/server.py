import socket, threading, os
import pydy_viz
class Server(threading.Thread):
    def __init__(self, port):
        threading.Thread.__init__(self)
        self.port = port
        host = ''
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.bind((host, self.port))
        self.socket.listen(5)

    def parse_data(self,data):
        static_path = os.path.dirname(pydy_viz.__file__)
        request = data.split(' ')[1][1:]
        file_path = os.path.join(static_path, request)
        try:
            send_buffer = open(file_path).read()
        except IOError:
            send_buffer = '404: file not found'
        
        return send_buffer        
        
    def listen_once(self):
        conn, addr = self.socket.accept()
        data = conn.recv(1024)
        sent_data = self.parse_data(data)
        conn.send(sent_data)
        return sent_data
        
        
    def run(self):
        host = ''
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.bind((host, self.port))
        print 'server started successfully, waiting...'
        s.listen(1)
        
        while 1:
            try:
                self.data = conn.recv(1024)
            except socket.error:
                print 'lost', addr, 'waiting..'
                s.listen(1)
                conn, addr = s.accept()
                print 'contact', addr, 'on', self.now()
                continue

            if not data:
                print 'lost', addr, 'waiting..'
                s.listen(1)
                conn, addr = s.accept()
            else:    
                sent_data = self.parse_data(self.data)
                conn.send(sent_data)


