"""
"""

from __future__ import absolute_import, division, print_function

import pytest
import numpy as np

from rxnnet import experiments



def test_constructor():
    expts = experiments.Experiments()
    assert list(expts.columns) == ['condition', 'varids', 'times']
    assert expts.nrow == 0


def test_add_experiment():
    expts = experiments.Experiments()
    expts.add_experiment(None, ['X'], [np.inf])
    assert expts.values.tolist() == [[None, ['X'], [np.inf]]]


def test_add():
    expts1 = experiments.Experiments()
    expts1.add_experiment(None, ['X'], [np.inf])

    expts2 = experiments.Experiments()
    expts2.add_experiment(('E','*',2), ['Y'], [1,2])

    expts = expts1 + expts2
    return expts

    assert expts.values.tolist() == [[None, ['X'], [np.inf]], 
                                     [('E','*',2), ['Y'], [1,2]]]
    assert expts.get_times() == [np.inf, 1.0, 2.0]
    assert expts.get_varids() == ['X', 'Y']
    assert expts.get_yids() == ['None, X, inf', 
                                '(E, *, 2), Y, 1', 
                                '(E, *, 2), Y, 2']




    