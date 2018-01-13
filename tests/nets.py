"""
"""

from __future__ import absolute_import, division, print_function
from collections import OrderedDict as OD

from rxnnet import network



def make_path2mah():
    """path2ma: A reaction network of two reactions forming a linear pathway, 
    where both reactions have mass-action-Haldane rate laws.
    """
    net = network.Network(id='path2mah')
    net.add_compartment(id='env')
    net.add_compartment(id='cell')
    net.add_species(id='C1', compartment='env', initial_value=2, 
                    is_constant=True)
    net.add_species(id='C2', compartment='env', initial_value=1, 
                    is_constant=True)
    net.add_species(id='X', compartment='cell', initial_value=0)
    net.add_reaction(id='R1', eqn='C1<->X', ratelaw='k1*(C1-X)', p={'k1':1})
    net.add_reaction(id='R2', eqn='X<->C2', ratelaw='k2*(X-C2)', p={'k2':2})
    net.compile()
    return net


def make_cycle3mah():
    """cycle3mah: A reaction network of three reactions forming a cycle,
    where all reactions have mass-action-Haldane rate laws.
    """
    net = network.Network(id='cycle3mah')
    net.add_compartment(id='env')
    net.add_compartment(id='cell') 
    net.add_species(id='C1', compartment='env', initial_value=2,
                    is_constant=True)  
    net.add_species(id='C2', compartment='env', initial_value=2,
                    is_constant=True)
    net.add_species(id='C3', compartment='env', initial_value=1,
                    is_constant=True)
    net.add_species(id='X1', compartment='cell', initial_value=1)
    net.add_species(id='X2', compartment='cell', initial_value=1)
    net.add_reaction(id='R1', eqn='C1+X1<->2 X2', ratelaw='k1*(C1*X1-X2**2)', 
                     p={'k1':1})
    net.add_reaction(id='R2', eqn='C2+X2<->X1', ratelaw='k2*(C2*X2-X1)', 
                     p={'k2':1})
    net.add_reaction(id='R3', eqn='X2<->C3', ratelaw='k3*(X2-C3)', 
                     p={'k3':1})
    net.compile()
    return net


def make_cycle4mah():
    """cycle4mah: A reaction network of four reactions forming a cycle, 
    where all reactions have mass-action-Haldane rate laws.
    """
    net = network.Network(id='cycle4mah')
    net.add_compartment(id='env')
    net.add_compartment(id='cell') 
    net.add_species(id='C1', compartment='env', initial_value=2,
                    is_constant=True)  
    net.add_species(id='C2', compartment='env', initial_value=2,
                    is_constant=True)
    net.add_species(id='C3', compartment='env', initial_value=1,
                    is_constant=True)
    net.add_species(id='X1', compartment='cell', initial_value=0)
    net.add_species(id='X2', compartment='cell', initial_value=0)
    net.add_species(id='X3', compartment='cell', initial_value=0)
    net.add_species(id='X4', compartment='cell', initial_value=0)
    net.add_reaction(id='R1', eqn='C1+X1<->2 X2')
    net.add_reaction(id='R2', eqn='X2+X3<->X1+X4')
    net.add_reaction(id='R3', eqn='C2+X4<->X3')
    net.add_reaction(id='R4', eqn='X2<->C3')
    return net


def make_path3mmh():
    """path3mmh: A reaction network of three reactions forming a linear pathway, 
    where all reactions have Michaelis-Menten-Haldane rate laws.

    seed=0; atol=rtol=1e-12; intermediate_output=True success =False fail
    """
    net = network.Network('path3mmh')
    net.add_compartment(id='env')
    net.add_compartment(id='cell')
    net.add_species(id='C1', compartment='env', initial_value=2, 
                         is_constant=True)
    net.add_species(id='C2', compartment='env', initial_value=0.5, 
                         is_constant=True)
    net.add_species(id='X1', compartment='cell', initial_value=1)
    net.add_species(id='X2', compartment='cell', initial_value=1)
    net.add_reaction(id='R1', eqn='C1<->X1', 
                          ratelaw='V1f/K1C1*(C1-X1/KE1)/(1+C1/K1C1+X1/K1X1)', 
                          p=OD([('V1f',1),('K1C1',1),('K1X1',1),('KE1',5)]))
    net.add_reaction(id='R2', eqn='X1<->X2', 
                          ratelaw='V2f/K2X1*(X1-X2/KE2)/(1+X1/K2X1+X2/K2X2)', 
                          p=OD([('V2f',1),('K2X1',1),('K2X2',1),('KE2',4)]))
    net.add_reaction(id='R3', eqn='X2<->C2', 
                          ratelaw='V3f/K3X2*(X2-C2/KE3)/(1+X2/K3X2+C2/K3C2)', 
                          p=OD([('V3f',1),('K3X2',1),('K3C2',1),('KE3',1)]))
    for KEid in ['KE1', 'KE2', 'KE3']:
        net.set_var_optimizable(KEid, False)
    net.compile()
    return net


path2mah = make_path2mah()
cycle3mah = make_cycle3mah()
cycle4mah = make_cycle4mah()
path3mmh = make_path3mmh()
