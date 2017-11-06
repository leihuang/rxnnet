# rxnnet

A Python package for representing **reaction networks**, and computing their **structures** and **behaviors**. Currently implemented as a wrapper of SloppyCell (http://sloppycell.sourceforge.net/). 

Mathematically, a reaction network can be represented as dx/dt = N v(x,p), where x, v and p are the concentration, rate and parameter vectors, respectively and N is the stoichiometry matrix.

Structures:
- stoichiometry matrix N
- left and right null spaces of N
- reduced stoichiometry matrix Nr (a selection of N rows with trivial left null space)
- link matrix L where N = L Nr

Behaviors:
- transients x(t)
- steady states s = x(inf) and J = v(s, p)
- parameter sensitivities of them dx(t)/dp, Rs = ds/dp and RJ = dJ/dp 
  - elasticities Ep = dv/dp and Ex = dv/dx 
  - control matrices Cs = -(N Es)^{-1} N and CJ = I + Es Cs

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
