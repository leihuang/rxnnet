"""
"""

from __future__ import division
import fractions
import re
import subprocess

import numpy as np
import pandas as pd

import util
reload(util)
from util import Matrix


def get_stoichiometry_matrix(net=None, rxnid2stoich=None, 
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
        integerize: if True, make all stoichcoefs integers; should be used with
                    caution as the ratelaws need to be changed accordingly to
                    keep the dynamics the same
    """
    if net:
        if only_dynvar:
            rowvarids = net.xids
        else:
            rowvarids = net.spids

        N = Matrix(np.zeros((len(rowvarids), net.vdim), dtype='int'),
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
            denoms = [fractions.Fraction(nz).limit_denominator().denominator 
                      for nz in nonzeros]
            denom = np.prod(list(set(denoms)))
            N.iloc[:,i] = col * denom
        
    return N


def get_steadystate_flux_matrix(net):
    """
    :param net:
    """
    N = net.N
    if hasattr(net, 'ss_flux_mat'):
        K = net.ss_flux_mat
        if np.allclose(N * K, [0]*net.xdim):
            return K
    elif N.rank == N.ncol:
        return Matrix(None, index=net.rxnids)
    else:
        # The following codes compute the INTEGER basis of right null space
        # of stoichiometry matrix.
        
        if np.allclose(net.N, net.N.astype('int')):
            # convert the matrix into a string recognizable by sage
            matstr = re.sub('\s|[a-z]|\(|\)', '', np.matrix(N).__repr__())
    
            # write a (sage) python script ".tmp_sage.py"
            # for more info of the sage commands: 
            # http://www.sagemath.org/doc/faq/faq-usage.html#how-do-i
            # -import-sage-into-a-python-script
            # http://www.sagemath.org/doc/tutorial/tour_linalg.html
            f = open('.tmp_sage.py', 'w')
            f.write('from sage.all import *\n\n')
            if np.allclose(net.N, net.N.astype('int')):
                f.write('A = matrix(ZZ, %s)\n\n' % matstr)  # integer field
            else:
                f.write('A = matrix(RR, %s)\n\n' % matstr)  # real field
            f.write('print kernel(A.transpose())')  # right nullspace vectors
            f.close()
            
            # call sage and run .tmp_sage.py
            out = subprocess.Popen(['sage', '-python', '.tmp_sage.py'],
                                   stdout=subprocess.PIPE)

            # process the output from sage
            vecstrs = out.communicate()[0].split('\n')[2:-1]
            #vecs = [eval(re.sub('(?<=\d)\s*(?=\d|-)', ',', vec)) 
            #        for vec in vecstrs]
            vecs = [filter(None, vec.strip('[]').split(' ')) for vec in vecstrs]
            try:
                vecs = [[int(elem) for elem in vec if elem] for vec in vecs]
            except ValueError:
                vecs = [[float(elem) for elem in vec if elem] for vec in vecs]
            fdids = ['J%d'%idx for idx in range(1, len(vecs)+1)]
            K = Matrix(np.transpose(vecs), net.rxnids, fdids)
            
        else:
            raise NotImplementedError

        net.ss_flux_mat = K

        return K


def get_pool_multiplicity_matrix(net):
    """Return a matrix whose row vectors are the multiplicities of 
    dynamic variables in the conservation pools. Mathematically, the matrix 
    has its rows spanning the left null space of the stoichiometry matrix.
    
    The function is computationally costly, because it calls *sage* to perform 
    matrix computations over the integer ring. 
    (Note that the matrix is converted to floats before being returned.)
    """
    N = net.N
    if hasattr(net, 'pool_mul_mat'):
        P = net.pool_mul_mat
        if np.allclose(P * N, [0]*net.vdim):
            return P
    elif N.rank == N.nrow:
        return Matrix(None, columns=net.xids)
    else:
        # The following codes compute the INTEGER basis of left null space
        # of stoichiometry matrix.

        # Convert the matrix into a string recognizable by sage.
        matstr = re.sub('\s|[a-z]|\(|\)', '', np.matrix(N).__repr__())
    
        ## Write a (sage) python script ".tmp_sage.py".
        # for more info of the sage commands: 
        # http://www.sagemath.org/doc/faq/faq-usage.html#how-do-i
        # -import-sage-into-a-python-script
        # http://www.sagemath.org/doc/tutorial/tour_linalg.html
        f = open('.tmp_sage.py', 'w')
        f.write('from sage.all import *\n\n')
        if np.allclose(net.N, net.N.astype('int')):
            f.write('A = matrix(ZZ, %s)\n\n' % matstr)  # integer field
        else:
            f.write('A = matrix(RR, %s)\n\n' % matstr)  # real field
        f.write('print A.kernel()')  # this returns the left nullspace vectors
        f.close()

        # Call sage and run .tmp_sage.py.
        out = subprocess.Popen(['sage', '-python', '.tmp_sage.py'],
                               stdout=subprocess.PIPE)
        
        # Process the output from sage.
        vecstrs = out.communicate()[0].split('\n')[2:-1]
        vecs = [eval(re.sub('(?<=\d)\s+(?=\d|-)', ',', vec)) for vec in vecstrs]

        poolids = ['Pool%d'%idx for idx in range(1, len(vecs)+1)]
        #poolids = net.xids[-len(vecs):][::-1]
        P = Matrix(vecs, poolids, N.rowvarids)

        # Clean things up: so far P can be, eg, 
        #        X1  X2  X3  X4
        # Pool1   0   0   1   1  # say, adenonine backbone
        # Pool2   2   1   3   2  # say, phospho group
        # We want it be the following, via Gaussian row reduction:
        #        X1  X2  X3  X4
        # Pool1   2   1   1   0
        # Pool2  -2  -1   0   1
        # so that X3 and X4 as the dependent dynvars can be easily 
        # selected and expressed as the linear combinations of 
        # independent dynvars

        # Adapted from here:
        # http://rosettacode.org/wiki/Reduced_row_echelon_form#Python
        def _to_reduced_row_echelon_form(m):
            M = m.values.tolist()
            lead = 0
            rowCount = len(M)
            columnCount = len(M[0])
            for r in range(rowCount):
                if lead >= columnCount:
                    return
                i = r
                while M[i][lead] == 0:
                    i += 1
                    if i == rowCount:
                        i = r
                        lead += 1
                        if columnCount == lead:
                            return
                M[i],M[r] = M[r],M[i]
                lv = M[r][lead]
                M[r] = [ mrx / float(lv) for mrx in M[r]]
                for i in range(rowCount):
                    if i != r:
                        lv = M[i][lead]
                        M[i] = [ iv - lv*rv for rv,iv in zip(M[r],M[i])]
                lead += 1
            return Matrix(M, m.rowvarids, m.colvarids)
        
        P = _to_reduced_row_echelon_form(P.ix[:, ::-1]).ix[::-1, ::-1]

    net.pool_mul_mat = P

    return P


def get_dependent_dynamical_variable_ids(net):
    """Return the dependent dynamical variables of the network (the dependence
    comes from conservation pools).

    For example, let P be:
            X1  X2  X3  X4
    Pool1   -2  -1   0   1
    Pool2    2   1   1   0
    
    It should return ['X3', 'X4']
    """
    P = net.P
    # dependent dynamic variables are picked at the end of each pool so that
    # networks that have been reordered or not will give the same variables
    if P is None:
        dxids = []    
    else:
        # pick the last xid in each pool
        dxids = [P.iloc[i][P.iloc[i]!=0].index[-1] for i in range(P.nrow)]
        
        # the following scenario should not be possible if P is in rref
        if len(dxids) != len(set(dxids)):
            # eg, ATP can be part of both the adenylate and phosphate pools;
            # in this case, it is easier to manually pick ddynvarids
            raise StandardError("The same dynamic variable has been chosen \
                from multiple pools as a dependent dynamic variable")
    return dxids


def get_independent_dynamical_variable_ids(net):
    """Return the independent dynamical variables of the network.
    """
    dxids = net.dxids
    return [xid for xid in net.xids if xid not in dxids]


def get_reduced_stoichiometry_matrix(net):
    """
    """
    ixids = net.ixids
    assert ixids == net.xids[:len(ixids)], "dynamical variables are not ordered"

    Nr = net.N.iloc[:len(ixids)]
    return Nr


def get_reduced_link_matrix(net):
    """
    L0: L = [I ]
            [L0]
    """
    ixids = net.ixids
    assert ixids == net.xids[:len(ixids)], "dynamical variables are not ordered"
    dxids = net.xids[len(ixids):]

    if len(dxids) == 0:
        L0 = Matrix(columns=net.xids)
    else:
        L0 = Matrix(-net.P.loc[:, ixids].values, dxids, ixids)
    return L0


def get_link_matrix(net):
    """
    L0: L = [I ]
            [L0]
            
    N = L * Nr
    """
    L0 = get_reduced_link_matrix(net)
    I = Matrix.eye(net.xids[:L0.ncol])
    L = pd.concat((I, L0))
    return L


def get_modules(net):
    raise NotImplementedError


