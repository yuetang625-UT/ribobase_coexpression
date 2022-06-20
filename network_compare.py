import networkx as nx
import netcomp as nc
import pandas as pd

def Graph_matrix(f,m):
    df=pd.read_csv("%s"%f,index_col=0)
    df=abs(df)
    df[df<m]=0
    df=df.replace(1,0)
    df=pd.DataFrame(df)
    df_np=df.to_numpy()
    G=nx.from_numpy_matrix(df_np)
    return G

G_RNA=Graph_matrix("rna_overlap_rho_clr.csv",0.5)
G_RIBO=Graph_matrix("ribo_overlap_rho_clr.csv",0.5)

A_RNA=nx.adjacency_matrix(G_RNA)
A_RIBO=nx.adjacency_matrix(G_RIBO)

deltacon=nc.deltacon0(A_RNA,A_RIBO)

#####https://github.com/peterewills/NetComp/blob/master/netcomp/linalg/fast_bp.py
from scipy import sparse as sps
import numpy as np
from numpy import linalg as la

def fast_bp(A,eps=None):
    """Return the fast belief propogation matrix of graph associated with A.
    Parameters
    ----------
    A : NumPy matrix or Scipy sparse matrix
        Adjacency matrix of a graph. If sparse, can be any format; CSC or CSR
        recommended.
    eps : float, optional (default=None)
        Small parameter used in calculation of matrix. If not provided, it is
        set to 1/(1+d_max) where d_max is the maximum degree.
    Returns
    -------
    S : NumPy matrix or Scipy sparse matrix
        The fast belief propogation matrix. If input is sparse, will be returned
        as (sparse) CSC matrix.
    Notes
    -----
    References
    ----------
    """
    n,m = A.shape
    ##
    ## TODO: implement checks on the adjacency matrix
    ##
    degs = np.array(A.sum(axis=1)).flatten()
    if eps is None:
        eps = 1/(1+max(degs))
    I = sps.identity(n)
    D = sps.dia_matrix((degs,[0]),shape=(n,n))
    # form inverse of S and invert (slow!)
    Sinv = I + eps**2*D - eps*A
    try:
        S = la.inv(Sinv)
    except:
        Sinv = sps.csc_matrix(Sinv)
        S = sps.linalg.inv(Sinv)
    return S

S_RNA_matrix = pd.DataFrame(S_RNA.toarray())