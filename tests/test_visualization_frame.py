from sympy.physics.mechanics import *
from sympy import symbols
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
inertia2 = inertia(B, Ixx[1], Iyy[1], Izz[1])

rigid_body1 = RigidBody('rigid_body1', P1, A, mass[0], (inertia1, P1))
rigid_body2 = RigidBody('rigid_body2', P2, B, mass[1], (inertia2, P2))

particle1 = Particle('particle1', P1, mass[0])                            

   
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
    
    
def test_vframe_with_rbody():
    
    frame2 = VisualizationFrame('frame2', rigid_body1, shape=mesh_shape1)
    
    assert frame2.name == 'frame2'
    assert frame1.reference_frame == A
    assert frame1.origin == P1
    assert frame1.shape == mesh_shape1
    
    frame1.name = 'frame2_'
    assert frame1.name == 'frame2_'
    
    frame1.reference_frame = B
    assert frame1.reference_frame == B
    
    frame1.origin = P2
    assert frame1.origin == P2
    
    frame1.shape = mesh_shape2
    assert frame1.shape == mesh_shape2    
    
    
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
