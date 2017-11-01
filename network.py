"""
"""

import re
from collections import OrderedDict as OD 

import SloppyCell.ReactionNetworks as scrn 


class Network(scrn.Network):
    """
    """
    def add_species(self, id, compartment, initial_value=0, *args, **kwargs):
        """A wrapper of SloppyCell.ReactionNetworks.Network.addSpecies.
        """
        self.addSpecies(id, compartment, initial_value, *args, **kwargs)
    
    
    def add_reaction(self, id, stoich=None, eqn=None, ratelaw=None, p=None, 
                     **kwargs):
        """A wrapper of SloppyCell.ReactionNetworks.Network.addReaction. 
        
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
        assert not (stoich is None) and (eqn is None)
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
        self.addReaction(id=rxnid, stoichiometry=stoich, kineticLaw=ratelaw, 
                         **kwargs)
    

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

