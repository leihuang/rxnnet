

from rxnnet import network


def test_make_path2mah():
    net = network.Network(id='net')
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

