"""
"""

import pandas as pd
import numpy as np


class Series(pd.Series):
    """
    """

    @property
    def _constructor(self):  
        return Series

        
    def randomize(self, seed=None, distribution='lognormal', **kwargs):
        """
        :param seed:
        :param distribution: name of the distribution used in randomization
            ('lognormal', 'normal', etc.). If 'lognormal': default is sigma=1
        :type distribution: str
        :param kwargs: sigma
        """
        if seed is not None:
            np.random.seed(seed)

        if distribution == 'lognormal':
            if 'sigma' not in kwargs:
                kwargs['sigma'] = 1
            x = self * np.random.lognormal(size=self.size, **kwargs)
        if distribution == 'normal':
            x = self + np.random.normal(size=self.size, **kwargs)
        return x 


class Matrix(pd.DataFrame):
    """
    """
    
    def __mul__(self, other):
        pass