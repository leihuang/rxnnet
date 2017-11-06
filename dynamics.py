"""
"""

import pandas as pd
from SloppyCell import daskr
from SloppyCell.ReactionNetworks import Dynamics


TOL = 1e-12

def integrate0(net, times, varids, p=None, tol=None, to_DF=False, use_daeint=False):
    if tol is None:
        tol = TOL
    
    if times[0] != 0:
        times = [0] + list(times)
        _prepend_zero = True
    else:
        _prepend_zero = False
        
    if p is not None:
        net.update_optimizable_vars(p)
   
    x0 = net.x0.copy()  # integration changes x0, hence the copying
    
    if not hasattr(net, 'res_function'):
        net.compile()
        
    if use_daeint:
        assert varids == net.xids
        out = daskr.daeint(res=net.res_function, t=times, y0=x0, yp0=[0]*net.xdim, 
                           atol=[tol]*net.xdim, rtol=[tol]*net.xdim, 
                           intermediate_output=False, rpar=net.constantVarValues)
        xt = out[0]
        xt[0] = net.x0  # somehow daskr.daeint messes up the first timepoint
    else:
        # Switch to Dynamics.integrate because of presumably easier indexing 
        # of requested varids.
        # Be careful of Dynamics.integrate's handling of times, esp. when 
        # the initial time is not zero.    
        traj = Dynamics.integrate(net, times, params=p, fill_traj=False,
                                  rtol=[tol]*net.xdim, atol=[tol]*net.xdim)
        xt = traj.copy_subset(varids).values
        
    if _prepend_zero:
        xt = xt[1:]
        times = times[1:]
    
    if to_DF:
        return butil.DF(xt, index=times, columns=varids)
    else:
        return xt


class Trajectory(pd.DataFrame):
    pass


def integrate(net, times, p=None, x0=None, tol=None, varids=None):
    """
    :param net:
    :type net:
    :param times:
    :type times: list or tuple
    :param varids: not implemented yet
    """
    if isinstance(times, tuple):
        intermediate_output=True
    else:
        intermediate_output=False

    t0 = times[0]
    if t0 != 0:
        _times = [0] + list(times)  # _times: for feeding the integrator
    else:
        _times = times

    if p is not None:
        net.p = p
    
    if x0 is None:
        x0 = net.x0

    if tol is None:
        tol = TOL
    
    if not hasattr(net, 'res_function'):
        net.compile()
        
    out = daskr.daeint(res=net.res_function, t=_times, 
                       y0=x0.copy(), yp0=[0]*net.xdim, 
                       atol=[tol]*net.xdim, rtol=[tol]*net.xdim, 
                       intermediate_output=intermediate_output, 
                       rpar=net.constantVarValues)

    traj = out[0]
    traj[0] = x0  # daskr.daeint messes up the first timepoint, or does it? FIXME
    times = out[1]

    if t0 != 0:
        idx_t0 = list(times).index(t0)  # used for dense output
        times = times[idx_t0:]
        traj = traj[idx_t0:]

    net.x = traj[-1]
    net.t = times[-1]

    return Trajectory(traj, index=pd.Index(times, name='time'), 
                      columns=net.xids)



