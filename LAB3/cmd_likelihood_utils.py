# ****************************************************************
#       Utility module for cmb_likelihood.py
# ****************************************************************
#
#   The module contains the following functions, usage being
#   returned_quantity = function(arguments)
#
#  *    N_cov    = get_noise_cov(rms)
#       F_cov    = get_foreground_cov(x,y,z)
#  *    C_ell    = get_C_ell_model(Q,n,lmax)
#       polys    = get_legendre_coeff(lmax)
#       P_ell_ij = get_legendre_mat(lmax,x,y,z)
#  *    S_cov    = get_signal_cov()
#  *    lnL      = get_lnL()
#
#
#  Routines marked by (*) are only partially implemented, and must be completed
#  before the program becomes functional.
#

import numpy as np
from scipy.special import legendre
import scipy.linalg as spl
import scipy.special as scis
import math

def get_noise_cov(sigma):
    return np.diag(sigma**2)


def get_foreground_cov(x,y,z):
    """
    Computing the foreground template covariance matrix, to marginalize over
    any monopole and dipole components in the map
    F_cov = large_value * sum(template_cov), where
    template_cov = np.outer(f, f^t). 
    For the monopole template, f is a constant.
    To account for a dipole of any orientation, we use each of the unit vector 
    components as a dipole template.
    """
    large_value = 1.0e3
    monopole = np.ones((len(x),len(x)))
    dipole = np.outer(x,x) + np.outer(y,y) + np.outer(z,z)
    return large_value * (monopole + dipole)

def get_C_ell_model(Q,n,lmax):
    """
    To be completed:
    Recursively compute a model power spectrum, C_ell, given the amplitude and
    spectral index parameters Q and n, on the range ell in [0,lmax],
    but with monopole and dipole terms set to 0.
    """
    # 1: Define array for power spectrum
    C_ell = np.zeros(lmax+1)
    # 2: Compute quadrupole (ell=2) term
    C_ell[2] = 4*np.pi/5*Q**2


    #Cl+1 = 4pi/5*q**2*(l- + n-1/2)gamma(l + n-1/2)gamma(9-n/2)/()
    # 3: Compute multipoles 3 through lmax recursively
    for l in range(3,lmax+1):
        C_ell[l] = (l - 1 + (n-1)/2.0)/(l - 1 + (5-n)/2.0)*C_ell[l-1]
        #C_ell[l] = C_ell[2]*math.gamma(l + (n-1)/2.0)*math.gamma((9-n)/2.0)/(math.gamma(l + (5-n)/2.0)*math.gamma((3+n)/2.0))

    return C_ell

def get_legendre_coeff(lmax):
    '''
    Helper routine for get_legendre_full. Computes Legendre polynomial
    coefficients for each multipole l, using scipy.special.legendre.
    Stores the result in a list of poly1d objects.
    Each such object returns the polynomial value when called with a 
    cos(theta) argument: P_ell = pol[l](costheta)
    '''
    leg = []
    for l in range(lmax+1):
        leg.append(legendre(l))
    return leg


def get_legendre_mat(lmax,x,y,z):
    '''
    Computing the full set of Legendre polynomial values needed to build the 
    signal covariance matrix.
    Uses helper function get_legendre_coeff for polynomial coefficients, and
    assembles a matrix of dimensions (ndata, ndata, lmax+1)
    '''
    leg = get_legendre_coeff(lmax)
    pos_vec = np.vstack([x,y,z]).T
    costheta =  np.dot(pos_vec,pos_vec.T)

    ndata = len(x)
    p_ell_ij = np.zeros((ndata,ndata,lmax+1))
    for l in range(lmax+1):
        p_ell_ij[:,:,l] = leg[l](costheta)
        
    return p_ell_ij


def get_signal_cov(C_ell, beam, pixwin, p_ell_ij):
    '''
    To be completed:
    Compute a (ndata,ndata) signal covariance matrix using the
    model power spectrum, instrument beam and pixel window function, and
    precomputed Legendre polynomials as input

    Hint: This can be done using a triple for-loop, but it is not necessary.
    Using NumPy array operations may get you a significant speed-up.
    '''
    lmax = len(p_ell_ij[0][0])-1
    
    vecell = 2.*np.arange(lmax+1, dtype=float) + 1.
    vecbeam = np.square(np.multiply(beam, pixwin))
    vec = np.multiply(vecell, vecbeam)
    vec = np.multiply(vec, C_ell)
    S_cov = np.dot(p_ell_ij, vec)
    """    
    lmax = len(p_ell_ij[0][0])-1
    
    B = np.zeros((len(p_ell_ij), len(p_ell_ij)))
    blpl = (beam*pixwin)**2

    for l in xrange(lmax+1):
        B += (2*l+1)*blpl[l]*C_ell[l]*p_ell_ij[:,:,l]
    """
    return S_cov/(4.*np.pi)


def get_lnL(data, cov):
    '''
    To be completed:
    Compute the quantity -2*lnL using the complete covariance matrix
    C = S+N+F, and the input data vector.

    Hint: This can be done directly by inverting the covariance matrix
    and computing a determinant, calling suitable NumPy or SciPy routines.
    Significant speedup can be gained by rather using the Cholesky decomposition
    method discussed in the project description notes.
    '''

    L = spl.cholesky(cov, lower = True)
    
    logdet = 2.*sum(np.log(np.diag(L)))

    x = spl.solve_triangular(L, data, lower=True)

    print logdet
    print np.dot(x.T,x)
    return np.dot(x.T,x) + logdet

