"""
"""

from __future__ import absolute_import, division, print_function

import numpy as np

from tests.nets import path2mah, cycle3mah, cycle4mah



def test_integrate():
    assert np.allclose(path2mah.integrate([0,1]), [[0.0], [1.2669505754952797]])
    assert np.allclose(path2mah.integrate([2,3]),
        [[1.330028330417476], [1.3331687869183433]])
    
    assert path2mah.integrate((0,1)).nrow > 10
    assert np.allclose(path2mah.integrate((0,1)).iloc[[0,-1],0], 
                       [0,1.2669505754952797])

    assert np.allclose(path2mah.integrate((1,2)).iloc[[0,-1],0],
                       [1.2669505754952797, 1.330028330418152])

