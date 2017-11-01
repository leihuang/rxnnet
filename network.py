"""
"""

import re
from collections import OrderedDict as OD 

import SloppyCell.ReactionNetworks as scrn 


class Network(scrn.Network):
    """
    """
    
    def add_reaction(self, id, stoich_or_eqn, ratelaw, p=None, **kwargs):
        """A convenient wrapper of SloppyCell.ReactionNetworks.Network.addReaction 
        
        Input:
            rxnid: a str; id of the reaction
            stoich_or_eqn: a mapping (stoich, eg, {'S':-1, 'P':1}) or 
                a str (eqn, eg, 'S<->P'); if an eqn is provided, 
                bidirectional arrow ('<->') denotes reversible reaction and 
                unidirectional arrow ('->') denotes irreversible reaction
            ratelaw:
            p: a dict; map from pid to pinfo, where pinfo can be 
                a float (pval) or a tuple (pval, is_optimizable); 
                eg, p={'V_R1':1, 'KM_S':2, 'KE_R1':(1, False)}
        """
        rxnid = id
        
        # get stoich
        if isinstance(stoich_or_eqn, str):
            eqn = stoich_or_eqn
            stoich = _eqn2stoich(eqn)
        else:
            stoich = stoich_or_eqn
        
        # add parameters
        if p is not None:
            for pid, pinfo in p.items():
                if isinstance(pinfo, tuple):
                    pval, is_optimizable = pinfo
                else:
                    pval, is_optimizable = pinfo, True
                if pid in self.parameters.keys():
                    if self.parameters.get(pid).value != pval:
                        raise ValueError("Value of parameter %s in reaction %s different."%\
                                         (pid, rxnid)) 
                else:
                    self.add_parameter(pid, pval, is_optimizable=is_optimizable)
    
        # add reaction
        self.addReaction(id=rxnid, stoichiometry=stoich, kineticLaw=ratelaw, 
                         **kwargs)
    

def _eqn2stoich(eqn):
    """Convert reaction equation (a str) to stoichiometry (a mapping).
    """
    def _unpack(s):
        # an example of s: ' 2 ATP '
        l = filter(None, s.split(' '))
        if len(l) == 1:
            sc_unsigned, spid = 1, l[0]   # sc: stoichcoef
        elif len(l) == 2:
            sc_unsigned, spid = eval(l[0]), l[1]
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

