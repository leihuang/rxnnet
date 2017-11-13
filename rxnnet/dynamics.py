"""
"""

from __future__ import division

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from SloppyCell import daskr


TOL = 1e-12


class Trajectory(pd.DataFrame):
    """
    """
    @property
    def _constructor(self):
        return Trajectory


    def __init__(self, *args, **kwargs):
        pd.DataFrame.__init__(self, *args, **kwargs)
        self.index.name = 'time'
        self.columns.name = 'varid'


    @property
    def times(self):
        return self.index.tolist()


    @property
    def varids(self):
        return self.columns.tolist()


    def plot(self, varids=None, figsize=None, filepath='', show=True):
        """
        """
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        if varids is not None:
            traj = traj[varids]
        else:
            traj = self
        ax.plot(traj.times, traj)
        ax.set_xlabel('time')
        ax.legend(traj.varids)
        ax.set_xlim(left=self.times[0])
        ax.set_ylim(bottom=0)
        if filepath:
            plt.savefig(filepath)
        if show:
            plt.show()



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

    t0 = times[0]
    if t0 == 0:
        _times = times
    else:
        _times = [0] + list(times)

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
    traj[0] = x0
    times = out[1]

    if t0 != 0:
        idx_t0 = list(times).index(t0)
        times = times[idx_t0:]
        traj = traj[idx_t0:]

    net.x = traj[-1]

    return Trajectory(traj, index=pd.Index(times, name='time'), 
                      columns=net.xids)



