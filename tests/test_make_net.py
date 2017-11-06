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

from rxnnet import network
reload(network)


@pytest.fixture()
def net():
    n = network.Network(id='net')
    n.add_compartment(id='env')
    n.add_compartment(id='cell')
    n.add_species(id='C1', compartment='env', initial_value=2, 
                  is_constant=True)
    n.add_species(id='C2', compartment='env', initial_value=1, 
                  is_constant=True)
    n.add_species(id='X', compartment='cell', initial_value=0)
    n.add_reaction(id='R1', eqn='C1<->X', ratelaw='k1*(C1-X)', p={'k1':1})
    n.add_reaction(id='R2', eqn='X<->C2', ratelaw='k2*(X-C2)', p={'k2':2})
    n.compile()
    return n


def test_get_N(net):
    assert net.N.values.tolist() == [[1.0, -1.0]]


def test_get_K(net):
    pass


def test_get_P(net):
    pass

    
def test_integrate(net):
    assert (net.integrate([0,1]).values.tolist() ==
        [[0.0], [1.2669505754952797]])
    assert (net.integrate([2,3]).values.tolist() ==
        [[1.330028330417476], [1.3331687869183433]])
    assert (net.integrate((0,1)).iloc[:3].values.tolist() ==
        [[0.0], [9.31322574614828e-13], [1.8626451492290054e-12]])
    assert (net.integrate((0,1)).iloc[-3:].values.tolist() ==
        [[1.2653340809755658], [1.2661904908822275], [1.2669505754952797]])
    assert (net.integrate((1,2)).iloc[:3].values.tolist() ==
        [[1.2669505754952797], [1.2670361148184137], [1.2678710886269715]])


def test_get_s():
    assert net.get_s().iloc[0] == 1.3333333333333333


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

