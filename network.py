"""
"""

import re
from collections import OrderedDict as OD 

import numpy as np
import SloppyCell.ReactionNetworks as scrn 

from rxnnet import util, dynamics
reload(dynamics)


class Network(scrn.Network, object):  # scrn.Network is an old-style class
    """
    """

    def add_species(self, id, compartment, initial_value=0, *args, **kwargs):
        """A wrapper of SloppyCell.ReactionNetworks.Network.add_species.
        """
        scrn.Network.add_species(self, id, compartment, initial_value, 
                                 *args, **kwargs)
    
    
    def add_reaction(self, id, stoich=None, eqn=None, ratelaw=None, p=None, 
                     **kwargs):
        """A wrapper of SloppyCell.ReactionNetworks.Network.add_reaction. 
        
        :param id: id of the reaction
        :type id: str
        :param stoich: a mapping (stoich, eg, {'S':-1, 'P':1}) or 
                a str (eqn, eg, 'S<->P'); if an eqn is provided, 
                bidirectional arrow ('<->') denotes reversible reaction and 
                unidirectional arrow ('->') denotes irreversible reaction
        :param eqn: 
        :type eqn: str
        :param ratelaw:
        :type ratelaw:
        :param p: map from pid to pinfo, where pinfo can be 
                a float (pval) or a tuple (pval, is_optimizable); 
                eg, p={'V_R1':1, 'KM_S':2, 'KE_R1':(1, False)}
        :type p:
        """
        rxnid = id
        
        # get stoich
        assert not (stoich is None and eqn is None)
        if eqn:
            stoich = _eqn2stoich(eqn)
            
        # add parameters
        if p is not None:
            for pid, pinfo in p.items():
                if isinstance(pinfo, tuple):
                    pval, is_optimizable = pinfo
                else:
                    pval, is_optimizable = pinfo, True
                self.add_parameter(pid, pval, is_optimizable=is_optimizable)
        
        if ratelaw is None:
            ratelaw = '0' 
        
        # add reaction
        scrn.Network.addReaction(self, id=rxnid, stoichiometry=stoich, 
                                 kineticLaw=ratelaw, **kwargs)


    @property
    def xdim(self):
        return len(self.dynamicVars)


    @property
    def pdim(self):
        return len(self.optimizableVars)


    @property
    def vdim(self):
        return len(self.reactions)

    
    @property
    def pids(self):
        return self.optimizableVars.keys()


    @property
    def p(self):
        return util.Series([var.value for var in self.optimizableVars], 
                           self.pids, dtype=np.float)

        
    @p.setter
    def p(self, p_new):
        self.set_var_vals(p_new.to_dict())


    @property
    def p0(self):
        return util.Series([var.initialValue for var in self.optimizableVars], 
                           self.pids, dtype=np.float)    


    @property
    def xids(self):
        return self.dynamicVars.keys()


    @property
    def x0(self):
        return util.Series([var.initialValue for var in self.dynamicVars], 
                           self.xids, dtype=np.float)


    @property
    def x(self):
        return util.Series([var.value for var in self.dynamicVars], 
                           self.xids, dtype=np.float)


    @x.setter
    def x(self, x_new):
        self.set_var_vals(x_new.to_dict())



    def integrate(self, *args, **kwargs):
        return dynamics.integrate(self, *args, **kwargs)
    integrate.__doc__ = dynamics.integrate.__doc__
    get_traj = integrate

    
    
    

def _eqn2stoich(eqn):
    """Convert reaction equation (a str) to stoichiometry (a mapping).
    """
    def _unpack(s):
        unpacked = filter(None, s.split(' '))  # an example of s: ' 2 ATP '
        if len(unpacked) == 1:
            sc_unsigned, spid = 1, unpacked[0]  # sc: stoichcoef
        elif len(unpacked) == 2:
            sc_unsigned, spid = eval(unpacked[0]), unpacked[1]
        else:
            raise 
        return spid, sc_unsigned
    
    # remove annotating species
    # eg, '(ST+)P->(ST+)G1P', where 'ST' (starch) is kept there to better
    # represent the chemistry
    eqn = re.sub('\(.*?\)', '', eqn)
    
    # re: '<?': 0 or 1 '<'; '[-|=]': '-' or '=' 
    subs, pros = re.split('<?[-|=]>', eqn)
    stoich = OD()
    
    if subs:
        for sub in subs.split('+'):
            subid, sc_unsigned = _unpack(sub)
            stoich[subid] = -1 * sc_unsigned
            
    if pros:
        for pro in pros.split('+'):
            proid, sc_unsigned = _unpack(pro)
            if proid in stoich:
                stoich[proid] = stoich[proid] + sc_unsigned
            else:
                stoich[proid] = sc_unsigned
        
    return stoich

