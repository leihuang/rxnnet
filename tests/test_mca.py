"""
"""

from __future__ import absolute_import, division, print_function

import numpy as np

from tests.nets import path2mah, cycle3mah, cycle4mah



cycle3mah.set_var_ics({'C1': 2, 'C2': 2, 'C3': 1})


def test_get_Ep():
    pass


def test_get_Es():
    assert np.allclose(cycle3mah.Es, [[2.0, -4.605551275463989], 
                                      [-1.0, 2.0], 
                                      [0.0, 1.0]])
    assert np.allclose(cycle3mah.nEs, [[5.070367516975986, -8.140735033951971],
                                       [-2.5351837584879977, 3.5351837584879977],
                                       [0.0, 1.7675918792439984]])

def test_get_jac():
    pass


def test_get_Cs():
    pass


def test_get_CJ():
    CJ_calc = [[0.2773500981126136, 0.5547001962252287, 0.16794970566215608],
               [0.2773500981126157, 0.5547001962252287, 0.16794970566215622],
               [0.27735009811261496, 0.5547001962252294, 0.16794970566215595]]
    assert np.allclose(cycle3mah.CJ, CJ_calc)
    assert np.allclose(cycle3mah.nCJ, CJ_calc)


def test_get_Rs():
    assert np.allclose(cycle3mah.Rs, 
        [[0.36132495094369227, 2.0254255396193805, -2.386750490563073],
         [0.3613249509436931, 0.7226499018873856, -1.0839748528310784]])


def test_get_RJ():
    pass