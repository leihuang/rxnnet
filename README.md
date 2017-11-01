# rxnnet

A Python package for representing **reaction networks** and simulating their behaviors (eg, transients, steady states and their parameter sensitivities). Currently implemented as a wrapper of SloppyCell (http://sloppycell.sourceforge.net/). 

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
