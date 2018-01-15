"""
"""

from __future__ import absolute_import, division, print_function
import itertools

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from SloppyCell import daskr
from SloppyCell.ReactionNetworks import Dynamics

from rxnnet import util
from infotopo import predict



ATOL = 1e-9
RTOL = 1e-9


class Trajectory(util.DF):
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


    def plot(self, varids=None, xylims=None, legendloc=None, 
             figsize=None, filepath='', show=True):
        """
        """
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        if varids is not None:
            traj = self[varids]
        else:
            traj = self
        ax.plot(traj.times, traj)
        ax.set_xlabel('time')
        ax.legend(traj.varids, loc=legendloc)
        
        if xylims is not None:
            if xylims[0] is not None:
                ax.set_xlim(xylims[0])
            if xylims[1] is not None:
                ax.set_ylim(xylims[1])
        else:
            ax.set_xlim(left=self.times[0])
            ax.set_ylim(bottom=0)

        if filepath:
            plt.savefig(filepath)
        if show:
            plt.show()



def integrate(net, times, p=None, x0=None, atol=None, rtol=None, varids=None):
    """A wrapper of SloppyCell.ReactionNetworks.Dynamics.integrate, where
    two previously nonintuitive behaviors are fixed here:
    1) x(t0) always falsely takes the value of x(t1), where t0 and t1 are 
        the first and second time points in the return traj;
    2) when the argument times does not start from 0, the integration goes
        from 0 to tmax-tmin.

    :param net:
    :param times: list or tuple
    :param atol:
    :param rtol:
    :param varids: 
    """
    if isinstance(times, tuple):
        intermediate_output = True
    else:
        intermediate_output = False

    t0 = times[0]
    if t0 == 0:
        _times = times
    else:
        _times = [0] + list(times)

    if p is not None:
        net.p = p

    if x0 is None:
        x0 = net.x0

    if atol is None:
        atol = ATOL
    if rtol is None:
        rtol = RTOL

    if not hasattr(net, 'res_function'):
        net.compile()

    """
    out = daskr.daeint(res=net.res_function, t=_times, 
                       y0=x0.copy(), yp0=[0]*net.xdim, 
                       atol=[atol]*net.xdim, rtol=[rtol]*net.xdim, 
                       intermediate_output=intermediate_output, 
                       rpar=net.constants)
    traj = out[0]
    traj[0] = x0
    times = out[1]
    """

    # Use SloppyCell.ReactionNetworks.Dynamics.integrate for now as it wraps 
    # around SloppyCell.daskr and is somehow more stable than daskr itself
    # but not much slower.
    # It automatically updates net.x hence no need to manually update x.
    out = Dynamics.integrate(net, _times, 
                             atol=[atol]*net.xdim, rtol=[rtol]*net.xdim, 
                             fill_traj=intermediate_output)

    if varids is not None:
        out = out.copy_subset(varids)

    traj = out.values
    times = out.timepoints
    varids = out.key_column.keys()

    if t0 != 0:
        idx_t0 = list(times).index(t0)
        times = times[idx_t0:]
        traj = traj[idx_t0:]

    return Trajectory(traj, index=pd.Index(times, name='time'), columns=varids)



def integrate_sensitivity(net, times, p=None, x0=None, rtol=None, 
                          varids=None):
    """

    :param net:
    :param times: list or tuple
    :param atol:
    :param rtol:
    :param varids: 
    """
    if isinstance(times, tuple):
        intermediate_output = True
    else:
        intermediate_output = False

    t0 = times[0]
    if t0 == 0:
        _times = times
    else:
        _times = [0] + list(times)

    if p is not None:
        net.p = p

    if x0 is None:
        x0 = net.x0

    if rtol is None:
        rtol = RTOL

    if not hasattr(net, 'res_function'):
        net.compile()

    out = Dynamics.integrate_sensitivity(net, _times, rtol=[rtol]*net.xdim, 
                                         fill_traj=intermediate_output)

    if varids is not None:
        if set(varids) < set(net.varids):
            varids = itertools.product(varids, net.pids)
        out = out.copy_subset(varids)

    traj = out.values
    times = out.timepoints
    varids = out.key_column.keys()

    if t0 != 0:
        idx_t0 = list(times).index(t0)
        times = times[idx_t0:]
        traj = traj[idx_t0:]

    traj = Trajectory(traj, index=pd.Index(times, name='time'), columns=varids)

    return traj



def get_predict(net, expts, **kwargs):
    """
    """
    # other cases have not been implemented yet
    assert expts.nrow == 1 and expts.get_conditions() == [None]

    varids = expts.get_varids()
    times = expts.get_times()

    atol = kwargs.pop('atol', None)
    rtol = kwargs.pop('rtol', None)
    
    def _f(p):
        traj = integrate(net, list(times), varids=varids, p=p, atol=atol, 
                         rtol=rtol)
        return traj.T.values.flatten()
    
    def _Df(p):  # note that sensitivity integration also returns x(t)
        straj = integrate_sensitivity(net, list(times), varids=varids, p=p, 
                                      rtol=rtol)
        return np.vstack([straj.values[:,i*net.pdim:(i+1)*net.pdim]
                          for i in range(len(varids))])

    pred = predict.Predict(f=_f, Df=_Df, p0=net.p0, pids=net.pids, 
                           yids=expts.get_yids(), expts=expts)
    
    return pred
