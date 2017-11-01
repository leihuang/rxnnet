

from rxnnet import network

def test_make_path2mah():
    net = network.Network(id='net')
    net.add_compartment(id='env')
    net.add_compartment(id='cell')
    net.add_species('C1', 'env', 2, is_constant=True)
    net.add_species('C2', 'env', 1, is_constant=True)
    net.add_species('X', 'cell', 0)
    net.add_reaction(id='R1', stoich_or_eqn='C1<->X', ratelaw='k1*(C1-X)', 
                     p={'k1':1})
    net.add_reaction(id='R2', stoich_or_eqn='X<->C2', ratelaw='k2*(X-C2)', 
                     p={'k2':2})
    net.compile()
    return net