from sympy.physics.mechanics import *
#from pydy-viz import VisualizationFrame

def test_visualization_frame():
	
	point_list1 = [[2, 3, 1], [4, 6, 2], \
                         [5, 3, 1], [5, 3, 6], \
                         [2, 8, 4], [7, 4, 1]]
    
    mesh_shape1 = MeshShape('mesh_shape1', points=point_list1, \
                               color='PINK', origin=[0,0,0])                     
    
    point_list2 = [[3, 1, 4], [3, 8, 2], \
                         [2, 1, 6], [2, 1, 1], \
                         [4, 1, 0], [1, 4, 3]]
    
    mesh_shape2 = MeshShape('mesh_shape2', points=point_list2, \
                               color='RED', origin=[0,0,0])                     
    
    #Testing constructor with rframe/point combination
    I = ReferenceFrame('I')
    O = Point('O')
    O.set_vel(I,0)
    frame1 = VisualizationFrame('frame1', [I, O], shape=mesh_shape1)
    
    assert frame1.name == 'frame1'
    assert frame1.reference_frame == I
    assert frame1.point == O
    assert frame1.shape == mesh_shape1
    
    frame1.name = 'frame1_'
    assert frame1.name == 'frame1_'
    
    #TODO : More tests for attributes
    
    
    
    #Testing constructor with rigidbody ..
    q, q1, Ixx, Iyy, Izz = dynamicsymbols('q q1 Ixx Iyy Izz')
    mass = symbols('mass')
    A = I.orientnew('A', 'Axia', [q, I.x])
    P = O.locatenew('P', q1 * A.y)
    P.v2pt_theory(O, I, A)
    inertia = inertia(A,Ixx,Iyy,Izz)
    rbody = RigidBody('rbody', P, A, mass, (inertia, O))
    
    frame2 = VisualizationFrame('frame2', rbody, shape=mesh_shape1)
    
    assert frame2.name == 'frame2'
    assert frame2.reference_frame == A
    assert frame2.point == P
    assert frame2.shape == mesh_shape1
    
    #TODO : more tests for properties
    
