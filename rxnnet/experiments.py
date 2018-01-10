"""
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import pandas as pd

from rxnnet import util

    

class Experiments(util.DF):
    """**Experiments** is a collection where each row corresponds to
    an individual **experiment**, represented by three objects:
        - condition: a 3-tuple of (target, type, strength) of **perturbation**;
            None for wildtype
        - varids:
        - times:
    """
    @property
    def _constructor(self):
        return Experiments
    
    
    def __init__(self, data=None, index=None, **kwargs):
        """
        """
        super(Experiments, self).__init__(data, index, 
            columns=['condition','varids','times'])
        self.index.name = 'experiment'

        for k, v in kwargs.items():
            setattr(self, k, v)
        

    def add(self, condition, varids, times):
        """
        :param condition: None or 3-tuple of (target, type, strength) of 
            perturbation
        :param varids: a list of variable ids
        :param times: a list of timepoints; finite values (dynamics) and inf
            (steady states) should not be in the same list
        """
        varids = [varids] if isinstance(varids, str) else list(varids)
        times = [times] if isinstance(times, float) else list(times)
        self.loc[self.nrow+1] = condition, varids, times


    def get_yids(self):
        """
        """
        def row2yids(row):
            condition, varids, times = row
            condid = '(%s)' %\
                ', '.join([''.join(map(str, tu)) for tu in condition])
            yids_row = [(condid, varid, time) for varid, time in 
                        util.get_product(varids, times)]
            return yids_row
        yids = [str(tu).replace("'", "")[1:-1] for tu in 
            util.flatten(self.apply(row2yids, axis=1), depth=1)]
        return yids


    def get_varids(self):
        """
        """
        return pd.Series(util.flatten(self.varids)).drop_duplicates().tolist()


    def get_times(self):
        """
        """
        return pd.Series(util.flatten(self.times)).drop_duplicates().tolist()

    
    def __add__(self, other):
        """
        """
        return Experiments(self.values.tolist()+other.values.tolist(),
                           index=range(1,self.nrow+other.nrow+1))


    def separate_by_time(self):
        """
        """
        assert all(self.times.apply(lambda times: np.inf not in times if 
            len(times)>1 else True)),\
            "times cannot have both finite values and inf"
        
        expts_dyn = expts[expts.times.apply(lambda times: times!=[np.inf])]
        expts_mca = expts[expts.times.apply(lambda times: times==[np.inf])]

        return expts_dyn, expts_mca



def get_experiments(varids, uids, us=None):
    """A convenience function for making experiments for reaction networks. 
    
    :param varids: 
    :param uids: id of control variables
    :param us: 
    """
    expts = Experiments()

    if uids == ['t']:
        expts.add(None, varids, us)
    else:
        conds = [tuple(zip(uids, ['=']*len(uids), u)) for u in us]
        for cond in conds:
            expts.add(cond, varids, [np.inf])
    return expts

