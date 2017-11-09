"""
"""

from __future__ import division

from SloppyCell import ExprManip as exprmanip

from util import Matrix


def get_Ex_str(net):
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


def get_Ep_str(net):
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


def get_param_elas_mat(net, p=None, normed=False):
    """
    TODO: compile or generate dynamic Python functions
    """
    net.update(p=p, t=np.inf)

    ns = net.namespace.copy()  # the namespace'd be contaminated without copy
    ns.update(net.varvals.to_dict())

    # assume no structural changes and attribute 'Ep_code' is up-to-date
    if not hasattr(net, 'Ep_code'):  
        Ep_code = get_Ex_str(net)[1]
    else:
        Ep_code = net.Ep_code

    Ep = Matrix(np.array(eval(Ep_code, ns)), net.rxnids, net.pids)

    if normed:
        nEp = Es.normalize(net.v, net.p)
        return nEp
    else:
        return Ep


def get_concn_elas_mat(net, p=None, normed=False):
    """
    """
    net.update(p=p, t=np.inf)

    ns = net.namespace.copy()  # the namespace'd be contaminated without copy
    ns.update(net.varvals.to_dict())

    # assume no structural changes and attribute 'Ex_code' is up-to-date
    if not hasattr(net, 'Ex_code'):  
        Ex_code = get_Ex_str(net)[1]
    else:
        Ex_code = net.Ex_code

    Es = Matrix(np.array(eval(Ex_code, ns)), net.rxnids, net.xids)

    if normed:
        nEs = Es.normalize(net.v, net.s)
        return nEs
    else:
        return Es


def get_jac_mat(net, p=None):
    """Return the jacobian matrix of the network at steady state. 

    *In the MCA context*, it is the jacobian of the independent vector field
    dxi/dt = Nr * v(xi,xd,p), so that the matrix is invertible.
    """
    net.update(p=p, t=np.inf)
    Nr, L, Es = net.Nr, net.L, net.Es
    M = Nr * Es * L
    return M


def get_concn_ctrl_mat(net, p=None, normed=False):
    """
    """
    net.update(p=p, t=np.inf)
    L, M, Nr = net.L, net.M, net.Nr
    Cs = -L * M.inv() * Nr
    if normed:
        nCs = Cs.normalize(net.s, net.J)
        return nCs
    else:
        return Cs
        

def get_flux_ctrl_mat(net, p=None, normed=False):
    """
    """
    net.update(p=p, t=np.inf)
    I, Es, Cs = Matrix.eye(net.rxnids, net.rxnids), net.Es, net.Cs
    CJ = I + Es * Cs
    if normed:
        nCJ = CJ.normalize(net.J, net.J)
        return nCJ
    else:
        return CJ
    

def get_concn_resp_mat(net, p=None, normed=False):
    """
    """
    net.update(p=p, t=np.inf)
    Ep, Cs = net.Ep, net.Cs
    Rs = Cs * Ep
    if normed:
        nRs = Rs.normalize(net.s, net.p)
        return nRs
    else:
        return Rs


def get_flux_resp_mat(net, p=None, normed=False):
    """
    """
    net.update(p=p, t=np.inf)
    Ep, CJ = net.Ep, net.CJ
    RJ = CJ * Ep
    if normed:
        nRJ = RJ.normalize(net.J, net.p)
        return nRJ
    else:
        return RJ

