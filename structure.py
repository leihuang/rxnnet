"""
"""

import fractions

import numpy as np

from util import Matrix


def get_stoich_mat(net=None, rxnid2stoich=None, 
                   only_dynvar=True, integerize=False):
    """Return the stoichiometry matrix (N) of the given network or 
    dict rxnid2stoich. Rows correspond to species, and columns correspond to 
    reactions.
    
    Input:
        rxnid2stoich: eg, {'R1':{'X1':1}, 'R2':{'X1':-1,'X2':1}}; 
                      net & rxnid2stoich: one and only one should be given
        only_dynvar: if True, use *dynamic* species as rows (keep out 
                        constant/buffered species);
                     if False, use species as row
        integerize: if True, make all stoichcoefs integers
    """
    if net:
        if only_dynvar:
            rowvarids = net.xids
        else:
            rowvarids = net.spids

        N = Matrix(np.zeros((len(rowvarids), net.vdim)),
                   rowvarids, net.rxnids)

        for spid in rowvarids:
            for rxnid in net.rxnids:
                try:
                    stoichcoef = net.rxns[rxnid].stoichiometry[spid]
                    if isinstance(stoichcoef, str):
                        stoichcoef = net.evaluate_expr(stoichcoef)
                    N.loc[spid, rxnid] = stoichcoef
                except KeyError:
                    pass  # mat[i,j] remains zero

    if rxnid2stoich:
        rxnids = rxnid2stoich.keys()
        spids = []
        for stoich in rxnid2stoich.values():
            for spid, stoichcoef in stoich.items():
                if int(stoichcoef) != 0 and spid not in spids:
                    spids.append(spid)
        N = Matrix(np.zeros((len(spids), len(rxnids))), spids, rxnids)
        for spid in spids:
            for rxnid in rxnids:
                try:
                    N.loc[spid, rxnid] = rxnid2stoich[rxnid][spid]
                except KeyError:
                    pass  # mat[i,j] remains zero
    
    # make all stoichcoefs integers by first expressing them in fractions
    if integerize: 
        for i in range(N.ncol):
            col = N.iloc[:,i]
            nonzeros = [num for num in butil.flatten(col) if num]
            denoms = [fractions.Fraction(nonzero).limit_denominator().denominator 
                      for nonzero in nonzeros]
            denom = np.prod(list(set(denoms)))
            N.iloc[:,i] = col * denom
        
    return N
