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
    @property
    def _constructor(self):
        return Matrix

    
    @property
    def rowvarids(self):
        return self.index.tolist()


    @property
    def colvarids(self):
        return self.columns.tolist()


    @property
    def nrow(self):
        return self.shape[0]


    @property
    def ncol(self):
        return self.shape[1]


    def __mul__(self, other):
        return Matrix(np.dot(self, other), self.index, other.columns)


    def inv(self):
        return Matrix(np.linalg.inv(self), self.columns, self.index)


    def normalize(self, y=None, x=None):
        """
        M = dy / dx
        M_normed = d logy/d logx = diag(1/y) * M * diag(x)
        """
        mat = self
        if y is not None:
            mat = Matrix.diag(1/y) * mat
        if x is not None:
            mat = mat * Matrix.diag(x)
        return mat 


    @property
    def rank(self):
        return np.linalg.matrix_rank(self)
        

    @staticmethod
    def eye(rowvarids, colvarids=None):
        if colvarids is None:
            colvarids = rowvarids
        return Matrix(np.eye(len(rowvarids)), rowvarids, colvarids)


    @staticmethod
    def diag(x):
        """Return the diagonal matrix of series x. 
        
        D(x) = diag(x)
        """
        return Matrix(np.diag(x), x.index, x.index)

