"""
    - is_ss
    - get_s_integration
    - get_s_rootfinding
    - get_s
    - set_ss
"""

from __future__ import division

import numpy as np
import scipy as sp

from util import Series


TMIN = 1e3
TMAX = 1e9
K = 1e3
TOL = 1e-10  # tolerance for testing steady state
METHOD = 'integration'


def get_dxdt(net):
    """Get the velocities of species dynamics, dx/dt. 
    (Another way of getting it is N*v.)
    """
    if not hasattr(net, 'res_function'):
        net.compile()
    dxdt = net.res_function(net.t, x, np.zeros(len(x)), net.constants)
    return Series(dxdt, index=net.xids)



def test_ss(net, tol=None):
    """Determine whether the net has reached steady state.
    
    Caution: The codes as they are implemented now would reach false
    conclusions in certain corner cases. For example, it will falsely conclude 
    that a net has reached steady state if the nets has a very long 
    characteristic time scale, which can happen when exploring some corners 
    of the parameter space, eg, k=1e-9. FIXME
    
    Output:
        True or False
    """
    if tol is None:
        tol = TOL
    if np.max(np.abs(net.get_dxdt())) < tol:
        return True
    else:
        return False
    

def get_s_integration(net, p=None, Tmin=None, Tmax=None, k=None, tol=None):
    """
    Input:
        T0:
        Tmax:
        k:
        tol:
    """
    # delayed argument binding to propagate changes in the settings of 
    # global variables
    if Tmin is None:
        Tmin = TMIN
    if Tmax is None:
        Tmax = TMAX
    if k is None:
        k = K
    if tol is None:
        tol = TOL
        
    if p is not None:
        net.update(p=p)
        
    if net.is_ss(tol=tol):
        return net.x
    else:
        t = Tmin
        while t <= Tmax:
            net.update(t=t)  # could spit out daeintException
            if net.is_ss(tol=tol_ss):
                return net.x
            else:
                t *= k
        raise Exception("Unable to reach steady state for p: %s"%\
                        str(net.p.tolist()))


def get_s_rootfinding(net):
    pass


def get_s(net):
    pass


def set_ss(net):
    pass

