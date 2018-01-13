"""
"""

from __future__ import absolute_import, division, print_function

import numpy as np

from tests.nets import path2mah, cycle3mah, cycle4mah, path3mmh



def test_dxdt():
    assert np.allclose(path3mmh.get_dxdt(), [0.2, 0.05])



def test_get_s():
    assert np.isclose(path2mah.get_s().iloc[0], 1.3333333333333333)
    assert np.allclose(cycle3mah.s, [3.302775637731996, 2.3027756377319952])




