"""
"""

from __future__ import division

import pandas as pd
import numpy as np
from SloppyCell import daskr


TOL = 1e-12


class Trajectory(pd.DataFrame):
    pass


def integrate(net, times, p=None, x0=None, tol=None, varids=None):
    """A wrapper of SloppyCell.daskr.daeint.

    Two nonintuitive behaviors of daskr.daeint are fixed here:
    1) x(t0) always falsely takes the value of x(t1), where t0 and t1 are 
        the first and second time points in the return traj;
    2) when the argument times does not start from 0, the integration goes
        from 0 to tmax-tmin.

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

    # t0 == 0 and x0 == None
    # t0 == 0 and x0 != None
    # t0 != 0 and x0 == None
    # t

    # t0 != 0, x0 provided



    # _times: for feeding the integrator
    t0 = times[0]
    if t0 == 0:
        x0 = net.x0
        _times = times
    elif t0 != 0 and x0 != net.x:
        x0 = net.x0
        _times = [0] + list(times)  
    else:  # t0 not 0 but t0 equal to net.t
        x0 = net.x
        _times = np.array(times) - t0

    if p is not None:
        net.p = p
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
    traj[0] = x0
    times = out[1]

    if t0 != 0 and t0 != net.t:
        idx_t0 = list(times).index(t0)
        times = times[idx_t0:]
        traj = traj[idx_t0:]
    if t0 != 0 and t0 == net.t:
        times = times + t0

    net.x = traj[-1]
    net.t = times[-1]  

    return Trajectory(traj, index=pd.Index(times, name='time'), 
                      columns=net.xids)



