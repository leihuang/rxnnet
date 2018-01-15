"""
"""

from __future__ import absolute_import, division, print_function
from collections import OrderedDict as OD 
import re
import copy

import numpy as np
import SloppyCell.ReactionNetworks as scrn
from SloppyCell import ExprManip as exprmanip

from rxnnet.util import Series
from rxnnet import structure, dynamics, steadystate, mca



class Network(object, scrn.Network):
    """

    scrn.Network is an old-style class and we want to make extensive use of 
    new-style class features such as properties, hence the double inheritance.
    """

    """
    __slots__ = ('add_species', 'add_reaction', 
                 'xdim', 'pdim', 'vdim',
                 'pids', 'p', 'p0', 
                 'xids', 'x', 'x0',
                 'c', 'constants', 
                 'rxnids', 'rxns', 'ratelaws',
                 'get_stoich_mat', 'N', 
                 'integrate', 'get_traj',
                 'get_s', 

                 't',

                 'set_var_vals',

        ) 
    """  


    def __init__(self, *args, **kwargs):
        scrn.Network.__init__(self, *args, **kwargs)

        # Thoughts on self.t:
        # Can easily introduce bugs, whenver the vector field is changed
        # but self.t and self.x0 are not reset
        # self.t = 0


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
        return Series([var.value for var in self.optimizableVars], 
                           self.pids, dtype=np.float)

        
    @p.setter
    def p(self, p_new):
        # It seems only calling SloppyCell's "update_optimizable_vars" (but 
        # not set_var_ic, set_var_vals, etc.) would udpate the network 
        # attribute "constantVarValues", which is important for many internal 
        # computations in SloppyCell.
        # I defined a new attribute "constants", which is implemented as a 
        # property and guarantees to stay up-to-date, with the intention of
        # replacing "constantVarValues", but still many internal computations
        # in SloppyCell are beyond reach.
        # Not handling this with enough care may introduce some cryptic bugs.
        self.update_optimizable_vars(p_new)


    @property
    def p0(self):
        return Series([var.initialValue for var in self.optimizableVars], 
                      self.pids, dtype=np.float)    


    @property
    def xids(self):
        return self.dynamicVars.keys()


    @property
    def spids(self):
        return self.species.keys()


    @property
    def Cids(self):
        return [spid for spid in self.spids if spid not in self.xids]


    @property
    def x0(self):
        return Series([var.initialValue for var in self.dynamicVars], 
                      self.xids, dtype=np.float)


    @property
    def x(self):
        return Series([var.value for var in self.dynamicVars], 
                      self.xids, dtype=np.float)


    @x.setter
    def x(self, x_new):
        try:
            self.set_var_vals(x_new.to_dict())
        except AttributeError:
            self.set_var_vals(dict(zip(self.xids, x_new)))


    @property
    def v(self):
        return Series(OD([(rxn.id, self.evaluate_expr(rxn.kineticLaw)) 
                          for rxn in self.reactions]), 
                      dtype=np.float)


    @property
    def c(self):
        return Series(OD([(sp.id, sp.value) for sp in self.species 
                          if sp.is_constant]), 
                      dtype=np.float)


    @property
    def constants(self):
        """Replacing the constantVarValues attribute of scrn.Network instances,
        as it is buggy and does not update upon changing p in certain ways.
        """
        return Series([var.value for var in self.constantVars], 
                      self.constantVars.keys(), dtype=np.float)


    @property
    def vars(self):
        return Series(OD(self.variables.items()), dtype=object)
    

    @property
    def varvals(self):
        return Series(OD([(var.id, var.value) for var in self.variables]), 
                      dtype=np.float)

    
    @property
    def varids(self):
        return self.variables.keys()


    @property
    def rxnids(self):
        return self.reactions.keys()


    @property
    def rxns(self):
        return Series(OD([(rxn.id, rxn) for rxn in self.reactions]), 
                      dtype=object)


    @property
    def ratelaws(self):
        return Series(OD([(rxn.id, rxn.kineticLaw) for rxn in self.reactions]),
                      dtype=object)

    
    @property
    def vids(self):
        return map(lambda rxnid: 'v_'+rxnid, self.rxnids)


    @property
    def Jids(self):
        return map(lambda rxnid: 'J_'+rxnid, self.rxnids)
    

    @property
    def asgrules(self):
        return Series(OD(self.assignmentRules.items()), dtype=object)


    def add_ratevars(self):
        """Add rate variables.
        """
        for rxn in self.reactions:
            rateid = 'v_' + rxn.id
            try:
                self.add_parameter(rateid, is_constant=False, 
                                   is_optimizable=False)
                self.add_assignment_rule(rateid, rxn.kineticLaw)
            except ValueError:
                pass


    def get_stoich_mat(self, *args, **kwargs):
        return structure.get_stoichiometry_matrix(net=self, *args, **kwargs)
    get_stoich_mat.__doc__ = structure.get_stoichiometry_matrix.__doc__


    @property
    def N(self):
        return structure.get_stoichiometry_matrix(self, only_dynvar=True, 
                                                  integerize=False)


    @property
    def K(self):
        return structure.get_steadystate_flux_matrix(self)


    @property
    def P(self):
        return structure.get_pool_multiplicity_matrix(self)


    @property
    def dxids(self):
        return structure.get_dependent_dynamical_variable_ids(self)


    @property
    def ixids(self):
        return structure.get_independent_dynamical_variable_ids(self)


    @property
    def Nr(self):
        return structure.get_reduced_stoichiometry_matrix(self)        


    @property
    def L0(self):
        return structure.get_reduced_link_matrix(self)


    @property
    def L(self):
        return structure.get_link_matrix(self)


    def integrate(self, *args, **kwargs):
        return dynamics.integrate(self, *args, **kwargs)
    integrate.__doc__ = dynamics.integrate.__doc__
    get_traj = integrate


    def get_dxdt(self):
        return steadystate.get_dxdt(self)
    get_dxdt.__doc__ = steadystate.get_dxdt.__doc__


    def test_ss(self, tol=None):
        return steadystate.test_ss(self, tol)
    test_ss.__doc__ = steadystate.test_ss.__doc__


    def set_ss(self, tol=None):
        steadystate.set_ss(self, tol)
    set_ss.__doc__ = steadystate.set_ss.__doc__


    def get_s(self, *args, **kwargs):
        return steadystate.get_s(self, *args, **kwargs)
    get_s.__doc__ = steadystate.get_s.__doc__


    def get_J(self, *args, **kwargs):
        return steadystate.get_J(self, *args, **kwargs)
    get_J.__doc__ = steadystate.get_J.__doc__


    @property
    def s(self):
        return self.get_s(method=steadystate.METHOD, tol=steadystate.TOL)


    @property
    def J(self):
        return self.get_J(method=steadystate.METHOD, tol=steadystate.TOL)


    @property
    def Es(self):
        return mca.get_concentration_elasticity_matrix(self)


    @property
    def nEs(self):
        return mca.get_concentration_elasticity_matrix(self, normed=True)


    @property
    def Ep(self):
        return mca.get_parameter_elasticity_matrix(self)


    @property
    def nEp(self):
        return mca.get_parameter_elasticity_matrix(self, normed=True)


    @property
    def M(self):
        return mca.get_jacobian_matrix(self)


    @property
    def Cs(self):
        return mca.get_concentration_control_matrix(self)


    @property
    def nCs(self):
        return mca.get_concentration_control_matrix(self, normed=True)

    @property
    def CJ(self):
        return mca.get_flux_control_matrix(self)


    @property
    def nCJ(self):
        return mca.get_flux_control_matrix(self, normed=True)


    @property
    def Rs(self):
        return mca.get_concentration_response_matrix(self)


    @property
    def nRs(self):
        return mca.get_concentration_response_matrix(self, normed=True)


    @property
    def RJ(self):
        return mca.get_flux_response_matrix(self)


    @property
    def nRJ(self):
        return mca.get_flux_response_matrix(self, normed=True)


    def update(self, p=None, x=None, t=None, **kwargs):
        """
        """
        if p is not None:
            self.p = p

        if x is not None:
            self.x = x

        if t is not None:
            if t == np.inf:
                tol = kwargs.pop('tol', steadystate.TOL)
                if not self.test_ss(tol=tol):
                    self.set_ss(tol=tol, **kwargs)
            else:
                raise NotImplementedError


    def replace_varid(self, varid_old, varid_new, only_expr=False):
        """Change id of reaction, species, or parameter.
        
        :param only_expr: bool; if True, only replace varid in expressions 
            such as reaction ratelaws or assignment rules but not variable ids; 
            useful when varid_new == 'varid_old * r'
        """
        if only_expr:
            f = lambda varid: varid
        else:
            f = lambda varid: varid_new if varid == varid_old else varid
        
        netid_new = f(self.id)
        net_new = self.__class__(netid_new)
        
        for var in self.variables:
            var_new = copy.deepcopy(var)
            var_new.id = f(var.id)
            net_new.variables.set(var_new.id, var_new)
        
        for rxn in self.reactions:
            rxn_new = copy.deepcopy(rxn)
            rxn_new.id = f(rxn.id)
            rxn_new.stoichiometry = OD(zip(map(f, rxn.stoichiometry.keys()),
                                           rxn.stoichiometry.values()))
            try:
                rxn_new.reactant_stoichiometry =\
                    OD(zip(map(f, rxn.reactant_stoichiometry.keys()),
                           rxn.reactant_stoichiometry.values()))
                rxn_new.product_stoichiometry =\
                    OD(zip(map(f, rxn.product_stoichiometry.keys()),
                           rxn.product_stoichiometry.values()))
            except AttributeError:  # some rxns have no reactant/product stoich
                pass

            rxn_new.parameters = set(map(f, rxn.parameters))
            rxn_new.kineticLaw = exprmanip.sub_for_var(rxn.kineticLaw, 
                                                       varid_old, varid_new)
            net_new.reactions.set(rxn_new.id, rxn_new)

        for varid, rule in self.assignmentRules.items():
            net_new.assignmentRules.set(f(varid), 
                exprmanip.sub_for_var(rule, varid_old, varid_new))
        
        for varid, rule in self.algebraicRules.items():
            net_new.algebraicRules.set(
                exprmanip.sub_for_var(varid, varid_old, varid_new), 
                exprmanip.sub_for_var(rule, varid_old, varid_new))

        for varid, rule in self.rateRules.items():
            net_new.rateRules.set(f(varid), 
                exprmanip.sub_for_var(rule, varid_old, varid_new))

        for eid, event in self.events.items():
            eid_new = f(eid)
            trigger_new = exprmanip.sub_for_var(event.trigger, varid_old, 
                                                varid_new)
            assignments_new = OD(zip(map(f, event.event_assignments.keys()),
                                     event.event_assignments.values()))
            net_new.add_event(eid_new, trigger_new, assignments_new)
            
        net_new.functionDefinitions = self.functionDefinitions.copy()
        for fid, f in net_new.functionDefinitions.items():
            fstr = 'lambda %s: %s' % (','.join(f.variables), f.math)
            net_new.namespace[fid] = eval(fstr, net_new.namespace)

        # _makeCrossReferences will take care of at least the following 
        # attributes:
        # assignedVars, constantVars, optimizableVars, dynamicVars, 
        # algebraicVars
        net_new._makeCrossReferences()
        return net_new


    def replace_varids(self, mapping):
        """
        :param mapping: a mapping from old varids to new varids 
        """
        net_new = self.copy()
        for varid_old, varid_new in mapping.items():
            net_new = net_new.replace_varid(varid_old, varid_new)
        return net_new


    def perturb(self, condition):
        """
        """
        if condition is None:  # wildtype  
            return self.copy()
        else:
            net = self.copy()
            if isinstance(condition[0], str):
                condition = (condition,)
            for perturbation in condition:
                varid, mode, strength = perturbation
                if mode in ['*', '/', '+', '-']:
                    varid_pert = '(%s%s%s)'%(varid, mode, strength)
                    net = net.replace_varid(varid, varid_pert, only_expr=True)  
                if mode == '=':
                    net.set_var_val(varid, strength)  # need to verify...
            return net


    def get_predict(self, expts, **kwargs):
        """
        """
        expts_dyn, expts_mca = expts.separate_by_time()

        if expts_dyn.nrow >= 1 and expts_mca.nrow == 0:
            return dynamics.get_predict(self, expts_dyn, **kwargs)
        elif expts_mca.nrow >= 1 and expts_dyn.nrow == 0:
            return mca.get_predict(self, expts_mca, **kwargs)
        else:
            return dynamics.get_predict(self, expts_dyn, **kwargs) +\
                mca.get_predict(self, expts_mca, **kwargs)


    def to_sbml(self, filepath):
        """
        """
        # Setting the following two attributes to None because otherwise 
        # _modifiers_ with stoichcoefs of 0 are not recognized 
        # when exporting the net to SBML, which would cause problems
        # when using the exported SBML in platforms such as Copasi or JWS
        for rxn in self.reactions:
            rxn.reactant_stoichiometry = None
            rxn.product_stoichiometry = None
        scrn.IO.to_SBML_file(self, filepath)



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

