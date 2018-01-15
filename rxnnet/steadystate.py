"""
    - is_ss
    - get_s_integration
    - get_s_rootfinding
    - get_s
    - set_ss
"""

from __future__ import absolute_import, division, print_function
import logging

import numpy as np
import scipy as sp

from rxnnet.util import Series, DF



TMIN = 1e3
TMAX = 1e9
K = 1e3
TOL = 1e-10  # tolerance for testing steady state
METHOD = 'integration'


logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


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
    dxdt = get_dxdt(net)

    logger.debug("x and dxdt for %s:\n%s" % (net.id, 
        DF([net.x, dxdt], ['x','dxdt']).T.__str__()))

    if np.max(np.abs(dxdt)) < tol:
        return True
    else:
        return False
    

def get_s_integration(net, p=None, Tmin=None, Tmax=None, k=None, 
                      tol=None, **kwargs):
    """
    :param T0:
    :param Tmax:
    :param k:
    :param tol:
    :param kwargs: kwargs for integrator
    """
    logger.debug("Get s through integration for network %s" % net.id)

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
            logger.debug("Integrate %s from %d to %d" %\
                (net.id, t_current, t_goal))

            _traj = net.integrate([0,t_goal-t_current], x0=x0, 
                                  **kwargs)

            if net.test_ss(tol=tol):
                return net.x
            else:
                t_current = t_goal
                t_goal *= k
                x0 = net.x

        logger.warn("Unable to reach steady state for p:\n%s"%\
                    net.p.__str__())
        raise Exception


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


def get_J(net, *args, **kwargs):
    """Get steady-state flux.
    """
    set_ss(net, *args, **kwargs)
    return Series(net.v.values, net.Jids)


def set_ss(net, tol=None, method=None, **kwargs):
    """
    """
    get_s(net, tol=tol, **kwargs)


