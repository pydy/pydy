from sympy.physics.mechanics import *
from sympy import symbols, sin, cos

#from pydy-viz import VisualizationFrame
#First we need to define some quantities... 
p = dynamicsymbols('p:3')
q = dynamicsymbols('q:3')

I = ReferenceFrame('I')
A = I.orientnew('A', 'body', p, 'XYZ') #p= [p1,p2,p3]
B = A.orientnew('B', 'body', q, 'XYZ')



O = Point('O')
P1 = O.locatenew('P1', 10*I.x + 10*I.y + 10*I.z)
P2 = P1.locatenew('P2', 10*I.x + 10*I.y + 10*I.z)


point_list1 = [[2, 3, 1], [4, 6, 2], [5, 3, 1], [5, 3, 6]]
point_list2 = [[3, 1, 4], [3, 8, 2], [2, 1, 6], [2, 1, 1]]

mesh_shape1 = MeshShape('mesh_shape1', points=point_list1, color='blue')    
mesh_shape2 = MeshShape('mesh_shape2', points=point_list2, color='red')     
                               
                               
Ixx = symbols('Ixx:2')
Iyy = symbols('Iyy:2')
Izz = symbols('Izz:2')
mass = symbols('mass:2')
    
inertia1 = inertia(A, Ixx[0], Iyy[0], Izz[0])

rigid_body = RigidBody('rigid_body1', P1, A, mass[0], (inertia1, P1))


particle = Particle('particle1', P1, mass[0])                            

#This is the transformation matrix from (A,P1) to (I,O)

transformation_matrix = [[cos(p[1])*cos(p[2]),  \
                sin(p[0])*sin(p[1])*cos(p[2]) + sin(p[2])*cos(p[0]), \
                sin(p[0])*sin(p[2]) - sin(p[1])*cos(p[0])*cos(p[2]), \
                10], [-sin(p[2])*cos(p[1]), \
                   -sin(p[0])*sin(p[1])*sin(p[2]) + \
                          cos(p[0])*cos(p[2]), \
                          sin(p[0])*cos(p[2]) + \
                  sin(p[1])*sin(p[2])*cos(p[0]),10], \
                   [sin(p[1]), -sin(p[0])*cos(p[1]), \
                                   cos(p[0])*cos(p[1]),10], \
                                    [0, 0, 0, 1]]
   
def test_vframe_with_rframe():
    frame1 = VisualizationFrame('frame1', [I, O], shape=mesh_shape1)
    
    assert frame1.name == 'frame1'
    assert frame1.reference_frame == I
    assert frame1.origin == O
    assert frame1.shape == mesh_shape1
    
    frame1.name = 'frame1_'
    assert frame1.name == 'frame1_'
    
    frame1.reference_frame = A
    assert frame1.reference_frame == A
    
    frame1.origin = P1
    assert frame1.origin == P1
    
    frame1.shape = mesh_shape2
    assert frame1.shape == mesh_shape2    
    
    assert frame1.transform(I, O).tolist() == transformation_matrix
    
def test_vframe_with_rbody():
    
    frame2 = VisualizationFrame('frame2', rigid_body, shape=mesh_shape1)
    
    assert frame2.name == 'frame2'
    assert frame2.reference_frame == A
    assert frame2.origin == P1
    assert frame2.shape == mesh_shape1
    
    frame2.name = 'frame2_'
    assert frame2.name == 'frame2_'
    
    frame2.reference_frame = B
    assert frame2.reference_frame == B
    
    frame2.origin = P2
    assert frame2.origin == P2
    
    frame2.shape = mesh_shape2
    assert frame2.shape == mesh_shape2    

    frame2.reference_frame = A
    frame2.origin = P1
    assert frame2.transform(I, O).tolist() == transformation_matrix
    
def test_vframe_with_particle():
    
    frame3 = VisualizationFrame('frame3', [particle1, A], shape=mesh_shape1)
    
    assert frame3.name == 'frame3'
    assert frame3.reference_frame == A
    assert frame3.origin == P1
    assert frame3.shape == mesh_shape1
    
    frame3.name = 'frame3_'
    assert frame3.name == 'frame3_'
    
    frame3.reference_frame = B
    assert frame3.reference_frame == B
    
    frame3.origin = P2
    assert frame3.origin == P2
    
    frame3.shape = mesh_shape2
    assert frame3.shape == mesh_shape2        

    frame3.reference_frame = A
    frame3.origin = P1
    assert frame3.transform(I, O).tolist() == transformation_matrix
