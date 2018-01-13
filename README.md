# rxnnet

A Python package for representing **reaction networks**, and computing their **structures** and **behaviors**. Currently the core of it is implemented as a wrapper of [SloppyCell](http://sloppycell.sourceforge.net/). 

Mathematically, a reaction network can be represented as dx/dt = N v(x,p), where x, v and p are the concentration, rate and parameter vectors, respectively and N is the stoichiometry matrix.

Structures:
- stoichiometry matrix N
- left and right null spaces of N
- reduced stoichiometry matrix N<sub>r</sub> (a selection of N's rows with trivial left null space) and link matrix L where N = L N<sub>r</sub>

Behaviors:
- dynamics x(t)
- steady states s = x(inf) and J = v(s, p)
- parameter sensitivities of them dx(t)/dp, R<sup>s</sup> = ds/dp and R<sup>J</sup> = dJ/dp 
  - elasticities E<sub>p</sub> = dv/dp and E<sub>x</sub> = dv/dx 
  - control matrices C<sup>s</sup> = -(N E<sub>s</sub>)<sup>-1</sup> N and C<sup>J</sup> = I + E<sub>s</sub> C<sup>s</sup>

<!---
Why a wrapper of SloppyCell: 
    - Coding styles
    - Extra functionalities such as steady states and mca

What rxnnet can do:
     - make networks
     - encode rate laws
     - network structures
     - simulate networks
     - get steady states
     - mca
-->

## Prerequisites
* [`SloppyCell`](http://sloppycell.sourceforge.net/)
* `pandas`
* `numpy`
* `matplotlib` (optional; for plotting dynamics)
* `scipy` (optional; for using rootfinding to get steady states)
* `sage` (optional; for [getting the integer-valued null spaces of a stoichiometry matrix](https://stackoverflow.com/questions/14407579/how-to-get-the-integer-eigenvectors-of-a-numpy-matrix))
* [`libsbml`](http://sbml.org/Software/libSBML/docs/python-api/) (optional; for importing and exporting SBML files)

## Usage examples

```python
import rxnnet

net = rxnnet.network.Network('net')
net.add_compartment(id='env')
net.add_compartment(id='cell')
net.add_species(id='C1', compartment='env', initial_value=2, is_constant=True)
net.add_species(id='C2', compartment='env', initial_value=1, is_constant=True)
net.add_species(id='X', compartment='cell', initial_value=0)
net.add_reaction(id='R1', eqn='C1<->X', ratelaw='k1*(C1-X)', p={'k1':1})
net.add_reaction(id='R2', eqn='X<->C2', ratelaw='k2*(X-C2)', p={'k2':2})

traj = net.integrate((0,10))
traj.plot()
```
