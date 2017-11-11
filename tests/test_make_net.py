"""
To test:

Nets:
    - path2mah
    - cyclemmh (using ratelaw module)

Behaviors:
    - p
    - traj
    - s (using both integration and rootfiding)
    - N, K
    - CJ, Cs, etc.
    - dxdt?
    
Design: 
    - Two setup functions via fixture?
    - Behavioral tests all in a class?
    - How to iterate over the (class) setup functions then?
"""

import pytest
import numpy as np

from rxnnet import network
reload(network)


@pytest.fixture()
def nets():
    net1 = network.Network(id='net1')
    net1.add_compartment(id='env')
    net1.add_compartment(id='cell')
    net1.add_species(id='C1', compartment='env', initial_value=2, 
                     is_constant=True)
    net1.add_species(id='C2', compartment='env', initial_value=1, 
                     is_constant=True)
    net1.add_species(id='X', compartment='cell', initial_value=0)
    net1.add_reaction(id='R1', eqn='C1<->X', ratelaw='k1*(C1-X)', p={'k1':1})
    net1.add_reaction(id='R2', eqn='X<->C2', ratelaw='k2*(X-C2)', p={'k2':2})
    net1.compile()

    net2 = network.Network(id='net2')
    net2.add_compartment(id='env')
    net2.add_compartment(id='cell') 
    net2.add_species(id='C1', compartment='env', initial_value=2,
                     is_constant=True)  
    net2.add_species(id='C2', compartment='env', initial_value=2,
                     is_constant=True)
    net2.add_species(id='C3', compartment='env', initial_value=1,
                     is_constant=True)
    net2.add_species(id='X1', compartment='cell', initial_value=0)
    net2.add_species(id='X2', compartment='cell', initial_value=0)
    net2.add_reaction(id='R1', eqn='C1+X1<->2 X2')
    net2.add_reaction(id='R2', eqn='C2+X2<->X1')
    net2.add_reaction(id='R3', eqn='X2<->C3')

    net3 = network.Network(id='net3')
    net3.add_compartment(id='env')
    net3.add_compartment(id='cell') 
    net3.add_species(id='C1', compartment='env', initial_value=2,
                     is_constant=True)  
    net3.add_species(id='C2', compartment='env', initial_value=2,
                     is_constant=True)
    net3.add_species(id='C3', compartment='env', initial_value=1,
                     is_constant=True)
    net3.add_species(id='X1', compartment='cell', initial_value=0)
    net3.add_species(id='X2', compartment='cell', initial_value=0)
    net3.add_species(id='X3', compartment='cell', initial_value=0)
    net3.add_species(id='X4', compartment='cell', initial_value=0)
    net3.add_reaction(id='R1', eqn='C1+X1<->2 X2')
    net3.add_reaction(id='R2', eqn='X2+X3<->X1+X4')
    net3.add_reaction(id='R3', eqn='C2+X4<->X3')
    net3.add_reaction(id='R4', eqn='X2<->C3')

    return [net1, net2, net3]


def test_get_N(nets):
    net1, net2, net3 = nets
    assert np.allclose(net1.N, [[1, -1]])
    assert np.allclose(net2.N, [[-1, 1, 0], [2, -1, -1]])


def test_get_P(nets):
    net1, net2, net3 = nets
    assert net1.P.shape == (0,1)
    assert net2.P.shape == (0,2)
    assert net3.P.shape == (1,4) and np.allclose(net3.P, [[0,0,1,1]])


def test_get_K(nets):
    net1, net2, net3 = nets
    assert net1.K.shape == (2,1) and np.allclose(net1.K, 1)
    assert net2.K.shape == (3,1) and np.allclose(net2.K, 1)


def test_get_Nr(nets):
    net1, net2, net3 = nets
    assert net3.ix    

    
def test_integrate(nets):
    net1, net2, net3 = nets
    assert all(np.isclose(net1.integrate([0,1]), [[0.0], [1.2669505754952797]]))
    assert all(np.isclose(net1.integrate([2,3]),
        [[1.330028330417476], [1.3331687869183433]]))
    assert all(np.isclose(net1.integrate((0,1)).iloc[:3],
        [[0.0], [9.31322574614828e-13], [1.8626451492290054e-12]]))
    assert all(np.isclose(net1.integrate((0,1)).iloc[-3:],
        [[1.2653340809755658], [1.2661904908822275], [1.2669505754952797]]))
    assert all(np.isclose(net1.integrate((1,2)).iloc[[0,-1]],
        [[1.2669505754952797], [1.330028330418152]]))


def test_get_s(nets):
    net1, net2, net3 = nets
    assert net1.get_s().iloc[0] == 1.3333333333333333


def test_get_Ep():
    pass


def test_get_Ex():
    pass


def test_get_jac():
    pass


def test_get_Cs():
    pass


def test_get_CJ():
    pass


def test_get_Rs():
    pass


def test_get_RJ():
    pass

