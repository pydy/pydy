from sympy.physics.mechanics import *
from sympy import sin, cos, symbols
#from pydy-viz import *

class TestVisualizationFrameScene(object):
    
    def __init__(self):
        #We define some quantities required for tests here..
        self.p = dynamicsymbols('p:3')
        self.q = dynamicsymbols('q:3')
        
        self.I = ReferenceFrame('I')
        self.A = self.I.orientnew('A', 'body', self.p, 'XYZ') #p= [p1,p2,p3]
        self.B = self.A.orientnew('B', 'body', self.q, 'XYZ')
        
        self.O = Point('O')
        self.P1 = self.O.locatenew('P1', 10 * self.I.x + \
                                      10 * self.I.y + 10 * self.I.z)
        self.P2 = self.P1.locatenew('P2', 10 * self.I.x + \
                                    10 * self.I.y + 10 * self.I.z)
        
        self.point_list1 = [[2, 3, 1], [4, 6, 2], [5, 3, 1], [5, 3, 6]]
        self.point_list2 = [[3, 1, 4], [3, 8, 2], [2, 1, 6], [2, 1, 1]]

        self.mesh_shape1 = MeshShape('mesh_shape1', \
                                points=self.point_list1, color='blue')    
        self.mesh_shape2 = MeshShape('mesh_shape2', \
                                points=self.point_list2, color='red')     
                               
                               
        self.Ixx, self.Iyy, self.Izz = symbols('Ixx Iyy Izz')
        
        self.mass = symbols('mass')
         
        self.inertia = inertia(self.A, self.Ixx, self.Iyy, self.Izz)

        self.rigid_body = RigidBody('rigid_body1', self.P1, self.A, \
                                 self.mass, (self.inertia, self.P1))
                          


        self.particle = Particle('particle1', self.P1, self.mass)                            

        self.transformation_matrix = [[cos(self.p[1])*cos(self.p[2]),  \
                  sin(self.p[0])*sin(self.p[1])*cos(self.p[2]) + \
                                     sin(self.p[2])*cos(self.p[0]), \
                  sin(self.p[0])*sin(self.p[2]) - \
                        sin(self.p[1])*cos(self.p[0])*cos(self.p[2]), \
                  10], [-sin(self.p[2])*cos(self.p[1]), \
                     -sin(self.p[0])*sin(self.p[1])*sin(self.p[2]) + \
                             cos(self.p[0])*cos(self.p[2]), \
                             sin(self.p[0])*cos(self.p[2]) + \
                     sin(self.p[1])*sin(self.p[2])*cos(self.p[0]),10], \
                    [sin(self.p[1]), -sin(self.p[0])*cos(self.p[1]), \
                                   cos(self.p[0])*cos(self.p[1]),10], \
                                     [0, 0, 0, 1]]
    

    def test_vframe_with_rframe(self):
        self.frame1 = VisualizationFrame('frame1', self.I, self.O, \
                                                shape=self.mesh_shape1)
    
        assert self.frame1.name == 'frame1'
        assert self.frame1.reference_frame == self.I
        assert self.frame1.origin == self.O
        assert self.frame1.shape is self.mesh_shape1
        
        self.frame1.name = 'frame1_'
        assert self.frame1.name == 'frame1_'
        
        self.frame1.reference_frame = self.A
        assert self.frame1.reference_frame == self.A
        
        self.frame1.origin = self.P1
        assert self.frame1.origin == self.P1
        
        self.frame1.shape = self.mesh_shape2
        assert self.frame1.shape is self.mesh_shape2    
        
        assert self.frame1.transform(self.I, self.O).tolist() == \
                                             self.transformation_matrix


    def test_vframe_with_rbody(self):
    
        self.frame2 = VisualizationFrame('frame2', self.rigid_body, \
                                                shape=self.mesh_shape1)
        
        assert self.frame2.name == 'frame2'
        assert self.frame2.reference_frame == self.A
        assert self.frame2.origin == self.P1
        assert self.frame2.shape == self.mesh_shape1
        
        self.frame2.name = 'frame2_'
        assert self.frame2.name == 'frame2_'
        
        self.frame2.reference_frame = self.B
        assert self.frame2.reference_frame == self.B
        
        self.frame2.origin = self.P2
        assert self.frame2.origin == self.P2
        
        self.frame2.shape = self.mesh_shape2
        assert self.frame2.shape is self.mesh_shape2    

        self.frame2.reference_frame = self.A
        self.frame2.origin = self.P1
        assert self.frame2.transform(self.I, self.O).tolist() == \
                                            self.transformation_matrix
                                            


    def test_vframe_with_particle(self):
        
        self.frame3 = VisualizationFrame('frame3', \
                                          self.particle1, self.A, \
                                                shape=self.mesh_shape1)
        
        assert self.frame3.name == 'frame3'
        assert self.frame3.reference_frame == self.A
        assert self.frame3.origin == self.P1
        assert self.frame3.shape is self.mesh_shape1
        
        self.frame3.name = 'frame3_'
        assert self.frame3.name == 'frame3_'
        
        self.frame3.reference_frame = self.B
        assert self.frame3.reference_frame == self.B
        
        self.frame3.origin = self.P2
        assert self.frame3.origin == self.P2
        
        self.frame3.shape = self.mesh_shape2
        assert self.frame3.shape is self.mesh_shape2        

        self.frame3.reference_frame = self.A
        self.frame3.origin = self.P1
        assert self.frame3.transform(self.I, self.O).tolist() == \
                                             self.transformation_matrix

    def test_vframe_without_name(self):
        self.frame4 = VisualizationFrame(self.I, self.O, \
                                               shape=self.mesh_shape1)
        
        assert self.frame4.name == 'unnamed'
        #To check if referenceframe and origin are defined 
        #properly without name arg
        assert self.frame4.reference_frame == self.I 
        assert self.frame4.origin == self.O
        assert self.frame4.shape is self.mesh_shape1
        
        self.frame4.name = 'frame1_'
        assert self.frame4.name == 'frame1_'
    
    def test_vframe_nesting(self):
        self.frame5 = VisualizationFrame('parent-frame'self.I, self.O, \
                                               shape=self.mesh_shape1)
        
        self.frame5.add_child_frames(frame1,frame2,frame3)                                                        
        assert self.frame5.get_children()[0] is self.frame1,
        assert self.frame5.get_children()[1] is self.frame2
        assert self.frame5.get_children()[2] is self.frame3
           
    def test_camera(self):
        #Camera is a subclass of VisualizationFrame, but without any
        #specific shape attached. We supply only ReferenceFrame,Point
        #to camera. and it inherits methods from VisualizationFrame
        camera = Camera('camera',I,O)
        
        assert camera.name == 'camera'
        assert camera.reference_frame == I
        assert cameta.origin == O
        
        camera.name = 'camera1'
        assert camera.name == 'camera1'
        
        camera.reference_frame = A
        assert camera.reference_frame == A
        
        camera.origin = P1
        assert camera.origin == P1
        
        #We can check transformation matrix for camera, in the extended 
        #example.
        
        #unnamed camera
        camera1 = Camera(I,O)
        assert camera1.name == 'unnamed'
        assert camera1.reference_frame == I
        assert cameta1.origin == O
