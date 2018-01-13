"""
To test:

Nets:
    - path2mah
    - cyclemmh (using ratelaw module)

Behaviors:
    - p
    - traj
    - s (using both integration and rootfiding)
    - N, K
    - CJ, Cs, etc.
    - dxdt?
    

"""

from __future__ import absolute_import, division, print_function

import numpy as np

from tests.nets import path2mah, cycle3mah, cycle4mah



def test_get_N():
    assert np.allclose(path2mah.N, [[1, -1]])
    assert np.allclose(cycle3mah.N, [[-1, 1, 0], [2, -1, -1]])


def test_get_P():
    assert path2mah.P.shape == (0,1)
    assert cycle3mah.P.shape == (0,2)
    assert cycle4mah.P.shape == (1,4) and np.allclose(cycle4mah.P, [[0,0,1,1]])


def test_get_K():
    assert path2mah.K.shape == (2,1) and np.allclose(path2mah.K, 1)
    assert cycle3mah.K.shape == (3,1) and np.allclose(cycle3mah.K, 1)


def test_get_Nr_L():
    assert np.allclose(cycle4mah.L * cycle4mah.Nr, cycle4mah.N)



