#*****************************************************************
# Chebyshev differentiation matrix
#  * Based on the formulation by Weideman and Reddy: 
#    https://web.iitd.ac.in/~mmehra/MATH4Q03/ft_gateway.cfm.pdf
#  * The following function is taken from `dmsuite`: 
#    https://github.com/ronojoy/pyddx/blob/master/sc/dmsuite.py), 
#    with small modifications. 
#-----------------------------------------------------------------

import numpy as np
import math as mt
import scipy.linalg as scla

def chebdif(N,M):
    """
    Chebyshev differentiation matrix up to max order M. 
    args:
       `M`: max order of differentiation
       `N`: number of nodes (in pyhysical space)
    return:   
       `D`: numpy array of shape (M,N,N)
    source: 
       https://github.com/ronojoy/pyddx/blob/master/sc/dmsuite.py
    """
    if M >= N:
        raise Exception('numer of nodes must be greater than M')        
    if M <= 0:
        raise Exception('derivative order must be at least 1')

    DM = np.zeros((M,N,N))    
    
    #n1 = (N/2); n2 = round(N/2.)     # indices used for flipping trick [Original]
    n1 = mt.floor(N/2); n2 = mt.ceil(N/2)     # indices used for flipping trick [Corrected]
    k = np.arange(N)                    # compute theta vector
    th = k*np.pi/(N-1)    

    # Compute the Chebyshev points
    #x = np.cos(np.pi*np.linspace(N-1,0,N)/(N-1))                # obvious way   
    x = np.sin(np.pi*((N-1)-2*np.linspace(N-1,0,N))/(2*(N-1)))   # W&R way
    x = x[::-1]
    
    # Assemble the differentiation matrices
    T = np.tile(th/2,(N,1))
    DX = 2*np.sin(T.T+T)*np.sin(T.T-T)               # trigonometric identity
    DX[n1:,:] = -np.flipud(np.fliplr(DX[0:n2,:]))    # flipping trick
    DX[range(N),range(N)]=1.                         # diagonals of D
    DX=DX.T

    C = scla.toeplitz((-1.)**k)           # matrix with entries c(k)/c(j)
    C[0,:]  *= 2
    C[-1,:] *= 2
    C[:,0] *= 0.5
    C[:,-1] *= 0.5

    Z = 1./DX                        # Z contains entries 1/(x(k)-x(j))
    Z[range(N),range(N)] = 0.        # with zeros on the diagonal.          

    D = np.eye(N)                    # D contains differentiation matrices.
    for ell in range(M):
        D = (ell+1)*Z*(C*np.tile(np.diag(D),(N,1)).T - D)      # off-diagonals    
        D[range(N),range(N)]= -np.sum(D,axis=1)        # negative sum trick
        DM[ell,:,:] = D                                # store current D in DM
    return x,DM