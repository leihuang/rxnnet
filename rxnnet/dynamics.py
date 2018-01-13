"""
"""

from __future__ import absolute_import, division, print_function

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from SloppyCell import daskr
from SloppyCell.ReactionNetworks import Dynamics

from rxnnet import util



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


    def plot(self, varids=None, figsize=None, filepath='', show=True):
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
        ax.legend(traj.varids)
        ax.set_xlim(left=self.times[0])
        ax.set_ylim(bottom=0)
        if filepath:
            plt.savefig(filepath)
        if show:
            plt.show()



def integrate(net, times, p=None, x0=None, atol=None, rtol=None, varids=None):
    """A wrapper of SloppyCell.daskr.daeint.

    Two nonintuitive behaviors of daskr.daeint are fixed here:
    1) x(t0) always falsely takes the value of x(t1), where t0 and t1 are 
        the first and second time points in the return traj;
    2) when the argument times does not start from 0, the integration goes
        from 0 to tmax-tmin.

    :param net:
    :param times: list or tuple
    :param atol:
    :param rtol:
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
    # We need to overwrite the Network attribute 'constantVarValues' 
    # in SloppyCell, as it somehow does not update upon changing parameters.
    net.constantVarValues = net.constants
    out = Dynamics.integrate(net, _times, 
                             atol=[atol]*net.xdim, rtol=[rtol]*net.xdim, 
                             fill_traj=intermediate_output)
    traj = out.values
    times = out.timepoints

    if t0 != 0:
        idx_t0 = list(times).index(t0)
        times = times[idx_t0:]
        traj = traj[idx_t0:]

    net.x = traj[-1]

    return Trajectory(traj, index=pd.Index(times, name='time'), 
                      columns=net.xids)


"""
def get_predict(net, expts, **kwargs):

    def get_xt(net, times, varids, p=None, tol=None, to_DF=False, use_daeint=False):
        if tol is None:
            tol = TOL
        
        if times[0] != 0:
            times = [0] + list(times)
            prepend_zero = True
        else:
            prepend_zero = False
            
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
            
        if prepend_zero:
            xt = xt[1:]
            times = times[1:]
        
        if to_DF:
            return util.DF(xt, index=times, columns=varids)
        else:
            return xt


    def get_dxtdp(net, times, varids, p=None, tol=None, to_DF=False):
        if tol is None:
            tol = TOL

        if times[0] != 0:
            times = [0] + list(times)
            prepend_zero = True
        else:
            prepend_zero = False
            
        straj0 = Dynamics.integrate_sensitivity(net, times, params=p, rtol=tol)
        straj = straj0.copy_subset(itertools.product(varids, net.pids))

        if prepend_zero:
            dat = np.vstack(np.hsplit(straj.values[1:], len(varids)))
            times = times[1:]
        else:
            dat = np.vstack(np.hsplit(straj.values, len(varids)))
            
        if to_DF:
            return util.DF(dat, index=itertools.product(varids, times), 
                            columns=net.pids)
        else:
            return dat
            

    assert expts.conds == [()], 'condition is not just wildtype.'
    #assert util.flatten(expts['varids']) == net.xids  # the restriction can be relaxed later
    varids = util.flatten(expts['varids'])
    
    
    def f(p):
        net.update(t=0)
        return get_xt(net, times util.flatten(expts['times']), p=p, varids=varids, **kwargs).T.flatten()
        
    
    def Df(p):
        net.update(t=0)
        return get_dxtdp(net, times util.flatten(expts['times']), p=p, varids=varids, **kwargs)

    pred = predict.Predict(f=f, Df=Df, p0=net.p0, pids=net.pids, 
                           yids=expts.yids, expts=expts)
    
    return pred
"""