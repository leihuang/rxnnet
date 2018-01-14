"""
"""

from __future__ import absolute_import, division, print_function

import numpy as np
from SloppyCell import ExprManip as exprmanip

from rxnnet import util
from infotopo import predict



def get_concentration_elasticity_string(net):
    """
    """
    Ex = []
    for rxnid in net.rxnids:
        ratelaw = exprmanip.sub_for_vars(net.rxns[rxnid].kineticLaw, 
                                         net.asgrules.to_dict())
        Ex_rxn = []
        for xid in net.xids:
            Ex_rxn.append(exprmanip.simplify_expr(
                exprmanip.diff_expr(ratelaw, xid)))
        Ex.append(Ex_rxn)
    Ex_str = str(Ex).replace("'", "")   
    Ex_code = compile(Ex_str, '', 'eval')  # compile to code object
    net.Ex_str, net.Ex_code = Ex_str, Ex_code
    return Ex_str, Ex_code



def get_parameter_elasticity_string(net):
    """
    """
    Ep = []
    for rxnid in net.rxnids:
        ratelaw = exprmanip.sub_for_vars(net.rxns[rxnid].kineticLaw, 
                                         net.asgrules.to_dict())
        Ep_rxn = []
        for pid in net.pids:
            Ep_rxn.append(exprmanip.simplify_expr(
                exprmanip.diff_expr(ratelaw, pid)))
        Ep.append(Ep_rxn)
    Ep_str = str(Ep).replace("'", "")   
    Ep_code = compile(Ep_str, '', 'eval')  # compile to code object
    net.Ep_str, net.Ep_code = Ep_str, Ep_code
    return Ep_str, Ep_code



def get_parameter_elasticity_matrix(net, p=None, normed=False, **kwargs):
    """
    TODO: compile or generate dynamic Python functions
    """
    net.update(p=p, t=np.inf, **kwargs)

    ns = net.namespace.copy()  # the namespace'd be contaminated without copy
    ns.update(net.varvals.to_dict())

    # assume no structural changes and attribute 'Ep_code' is up-to-date
    if not hasattr(net, 'Ep_code'):  
        Ep_code = get_parameter_elasticity_string(net)[1]
    else:
        Ep_code = net.Ep_code

    Ep = util.Matrix(np.array(eval(Ep_code, ns)), net.rxnids, net.pids)

    if normed:
        nEp = Es.normalize(net.v, net.p)
        return nEp
    else:
        return Ep
get_Ep = get_parameter_elasticity_matrix



def get_concentration_elasticity_matrix(net, p=None, normed=False, **kwargs):
    """
    """
    net.update(p=p, t=np.inf, **kwargs)

    ns = net.namespace.copy()  # the namespace'd be contaminated without copy
    ns.update(net.varvals.to_dict())

    # assume no structural changes and attribute 'Ex_code' is up-to-date
    if not hasattr(net, 'Ex_code'):  
        Ex_code = get_concentration_elasticity_string(net)[1]
    else:
        Ex_code = net.Ex_code

    Es = util.Matrix(np.array(eval(Ex_code, ns)), net.rxnids, net.xids)

    if normed:
        nEs = Es.normalize(net.v, net.s)
        return nEs
    else:
        return Es
get_Es = get_concentration_elasticity_matrix



def get_jacobian_matrix(net, p=None, **kwargs):
    """Return the jacobian matrix of the network at steady state. 

    *In the MCA context*, it is the jacobian of the independent vector field
    dxi/dt = Nr * v(xi,xd,p), so that the matrix is invertible.
    """
    net.update(p=p, t=np.inf, **kwargs)
    Nr, L, Es = net.Nr, net.L, net.Es
    M = Nr * Es * L
    return M



def get_concentration_control_matrix(net, p=None, normed=False, **kwargs):
    """
    """
    net.update(p=p, t=np.inf, **kwargs)
    L, M, Nr = net.L, net.M, net.Nr
    Cs = -L * M.inv() * Nr
    if normed:
        nCs = Cs.normalize(net.s, net.J)
        return nCs
    else:
        return Cs
get_Cs = get_concentration_control_matrix

        

def get_flux_control_matrix(net, p=None, normed=False, **kwargs):
    """
    """
    net.update(p=p, t=np.inf, **kwargs)
    I, Es, Cs = util.Matrix.eye(net.rxnids, net.rxnids), net.Es, net.Cs
    CJ = I + Es * Cs
    if normed:
        nCJ = CJ.normalize(net.J, net.J)
        return nCJ
    else:
        return CJ
get_CJ = get_flux_control_matrix

    

def get_concentration_response_matrix(net, p=None, normed=False, **kwargs):
    """
    """
    net.update(p=p, t=np.inf, **kwargs)
    Ep, Cs = net.Ep, net.Cs
    Rs = Cs * Ep
    if normed:
        nRs = Rs.normalize(net.s, net.p)
        return nRs
    else:
        return Rs
get_Rs = get_concentration_response_matrix



def get_flux_response_matrix(net, p=None, normed=False, **kwargs):
    """
    """
    net.update(p=p, t=np.inf, **kwargs)
    Ep, CJ = net.Ep, net.CJ
    RJ = CJ * Ep
    if normed:
        nRJ = RJ.normalize(net.J, net.p)
        return nRJ
    else:
        return RJ
get_RJ = get_flux_response_matrix



def get_predict(net, expts, **kwargs):
    """
    :param net:
    :param expts:
    :param kwargs: kwargs for networks getting steady states
    """
    assert expts.get_times() == [np.inf]

    # require varids to be the same; have not implemented heterogeneous varids 
    # yet, but it can otherwise be easily achieved through predict composition
    assert all(expts.varids.apply(lambda varids: varids==expts.varids.iloc[0]))

    varids = expts.get_varids()

    if set(varids) <= set(net.xids):
        vartype = 's'
        idxs = [net.xids.index(varid) for varid in varids]
    elif set(varids) <= set(net.Jids):
        vartype = 'J'
        idxs = [net.Jids.index(varid) for varid in varids]
    else:
        vartype = 'sJ'
        idxs = [(net.xids+net.Jids).index(varid) for varid in varids]
        
    net = net.copy()
    if not net.compiled:
        net.compile()

    nets = [net.perturb(cond) for cond in expts.condition]
    for net in nets:
        _Ex_str = get_concentration_elasticity_string(net)
        _Ep_str = get_parameter_elasticity_string(net)
        _L = net.L
        _Nr = net.Nr

    def _f(p):
        y = []
        for net in nets:
            net.update(p=p, t=np.inf, **kwargs)
            if vartype == 's':
                y_cond = net.x[idxs]
            if vartype == 'J':
                y_cond = net.v[idxs]
            if vartype == 'sJ':
                y_cond = np.concatenate((net.x, net.v))[idxs]
            y.extend(y_cond.tolist())
        return np.array(y)
    
    def _Df(p):
        jac = []
        for net in nets:
            net.update(p=p, t=np.inf, **kwargs)
            if vartype == 's':
                #jac_cond = get_Rs(net, p, Nr, L, to_mat=1, **kwargs).loc[varids].dropna()  # why dropna?
                jac_cond = net.Rs.iloc[idxs]
            if vartype == 'J':
                jac_cond = net.RJ.iloc[idxs]
            if vartype == 'sJ': 
                jac_cond = np.vstack((net.Rs, net.RJ)).iloc[idxs]
            jac.extend(jac_cond.values.tolist())
        return np.array(jac)

    
    pred = predict.Predict(f=_f, Df=_Df, p0=net.p0, pids=net.pids, 
                           yids=expts.get_yids(), expts=expts, nets=nets)
    
    return pred