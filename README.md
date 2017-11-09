# rxnnet

A Python package for representing **reaction networks**, and computing their **structures** and **behaviors**. Currently implemented as a wrapper of [SloppyCell](http://sloppycell.sourceforge.net/). 

Mathematically, a reaction network can be represented as dx/dt = N v(x,p), where x, v and p are the concentration, rate and parameter vectors, respectively and N is the stoichiometry matrix.

Structures:
- stoichiometry matrix N
- left and right null spaces of N
- reduced stoichiometry matrix N<sub>r</sub> (a selection of N rows with trivial left null space)
- link matrix L where N = L N<sub>r</sub>

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
* `scipy` (if using rootfinding to get steady states)
* `sage` (for [getting the null spaces of a stoichiometry matrix](https://stackoverflow.com/questions/14407579/how-to-get-the-integer-eigenvectors-of-a-numpy-matrix))

## Usage examples

```python
import rxnnet

net = rxnnet.Network()
net.add_species('C1')
net.add_species('X')
net.add_species('C2')
net.add_reaction('C1<->X')
net.add_reaction('X<->C2')

traj = net.get_traj((0,10))
traj.plot()
```
