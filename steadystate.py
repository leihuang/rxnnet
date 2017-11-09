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


TMIN = 1e2
TMAX = 1e8
K = 1e2
TOL = 1e-10  # tolerance for testing steady state
METHOD = 'integration'


def get_dxdt(net):
    """Get the velocities of species dynamics, dx/dt. 
    (Another way of getting it is N*v.)
    """
    if not hasattr(net, 'res_function'):
        net.compile()
    # assume it is an autonomous dynamical system 
    # (time-invariant, hence t=0 here)
    dxdt = net.res_function(0, net.x, np.zeros(net.xdim), net.constants)
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
    if np.max(np.abs(get_dxdt(net))) < tol:
        return True
    else:
        return False
    

def get_s_integration(net, p=None, Tmin=None, Tmax=None, k=None, 
                      tol=None, tol_integrator=None):
    """
    :param T0:
    :param Tmax:
    :param k:
    :param tol:
    :param tol_integrator:
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
        net.p = p
        
    if net.test_ss(tol=tol):
        return net.x
    else:
        t_current= 0
        t_goal = Tmin
        x0 = net.x0

        while t_goal <= Tmax:
            _traj = net.integrate([0,t_goal-t_current], x0=x0, 
                                  tol=tol_integrator)

            if net.test_ss(tol=tol):
                return net.x
            else:
                t_goal *= k
                t_current = t_goal
                x0 = net.x

        raise Exception("Unable to reach steady state for p: %s"%\
                        str(net.p.tolist()))


def get_s_rootfinding(net):
    pass


def get_s(net, method=None, *args, **kwargs):
    """
    """
    if method is None:
        method = METHOD

    if method == 'integration':
        return get_s_integration(net, *args, **kwargs)
    elif method == 'rootfinding':
        return get_s_rootfinding(net, *args, **kwargs)
    else:
        raise ValueError


def set_ss(net, tol=None, **kwargs):
    """
    """
    get_s(net, tol=tol, **kwargs)


