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


def test_get_traj(net):
    assert net.id == 'net'


def test_get_s():
    pass


def test_get_N():
    pass


def test_get_K():
    pass
