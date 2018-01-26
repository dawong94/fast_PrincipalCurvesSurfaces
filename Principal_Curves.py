"""Principal_Curves_Surfaces"""

# Author: James McQueen <jmcq@u.washington.edu>
#         Xiao Wang <wang19@u.washington.edu>
# LICENSE: Simplified BSD https://github.com/mmp2/megaman/blob/master/LICEN
import numpy as np
import scipy as sci
class Principal_Curves():
    def __init__(self, d, h,tol):
        self.d = d
        self.h = h
        self.tol = tol
                
    def fit(data, points):
        self.points_ = Principal_Curve(data, points, self.tol,self.d,self.h)
        
    def fit_transform(data, points):
        self.fit(data, points)
        return self.points_    
        
def constant(h, n, N):
    const = (1.0/N)*((h)**(-n))*(2.0*np.pi)**(-n/2.0)
    return const
            
        
def principal_curve_surface(data, points, tol, d, h):
    """    
    Principal curves are generally defined as smooth curves passing
    through the middle of the data point. This function will return
    every point on the principal curve or surface which is local maximum 
    of the probablity density in the local orthogonal subspace.
    
    Parameters    
    --------------------
    data: array-like, shape= (n_samples,N_samples)
    points: array-like, shape= (m_samples,M_samples)
    tol: float 
        tolerance   
    d: integer
        demension of principal surface
    h: float
        kernel bandwidth
    
    Returns
    -------
    points: array-like, shape (n_samples,N_samples)    
            projection points
    """
    data = np.atleast_2d(data)
    n, N = data.shape
    points = np.atleast_2d(points)
    m, M = points.shape
    if m == 1 and M == n: # row vector
        points = np.reshape(points, (n, 1))
        m, M = points.shape
    for k in range(M):
        
        points[:, k]= projectpoint(data, points[:,k], tol, d, h)   
#    print "xiao",points 
    return points
    

def projectpoint(data, point, tol, d, h):
    "Calculating the projection points"
    n, N = data.shape
    converged = False
    iter = 0
    while not(converged):        
        proj = get_projected_point(data, point, n, d, h)
        diff = np.linalg.norm(np.reshape(point, (n, 1)) - proj)
        point = np.reshape(proj, (n, ))
        iter = iter +1
        if diff < tol:
            converged = True
        if iter > 500:
            converged = True
            print "maximum iterations exceeded"
    return point 
    
def get_projected_point(data, point, n, d, h):
    n, N = data.shape
    
    H = hessian(data, point,h)   
    const = constant(h, n, N)
    w, v = np.linalg.eigh(H)  
    
    index = np.argsort(w) # arguments that sort from small to large
    V = np.reshape(v[:, index[range(n-d)]], (n, d))
    
    
#    print "eig",V
#"================================================================================"    
#    n, N = data.shape
#    c,u = C(data, point, h) 
#    const = constant(h, n, N)
#    U1 =  u*np.sqrt(c)
#    U,s,b= np.linalg.svd(U1/np.sqrt(N), full_matrices=False)
#    
#    s =  (s**2  - (np.sum(c)/(h**2))/N)
# #   print "svd",U,s
#    
#    index = np.argsort(s) # arguments that sort from small to large
#    V = np.reshape(U[:, index[range(n-d)]], (n, d))
##    print "svd", V    
#"================================================================================="    
    
    ospace = np.dot(V, V.T)
   
    proj = np.reshape(mean_shift(data, point, h), (n, 1)) - np.reshape(point, (n, 1))
    proj = np.dot(ospace, proj) + np.reshape(point, (n, 1))
    #proj = proj + np.reshape(point, (n, 1))
    return proj
    
def mean_shift(data, x, h):
    "Calculate the mean-shift"
    data = np.atleast_2d(data)
    n, N = data.shape
    const = constant(h, n, N)
    x = np.atleast_2d(x)
    m, M = x.shape
#    if M != 1:
#        x = np.reshape(x, (n, 1))
    c, u = C(data, x, h)    
    m_x = np.sum(c*data, axis = 1)/ np.sum(c)  
  
    m_x = np.reshape(m_x, (n, 1))
    return m_x
    


def hessian(data, x, h):
    "Calculate the Hessian"
    data = np.atleast_2d(data)
    n, N = data.shape   
    x = np.atleast_2d(x)
    m, M = x.shape
    if M != 1:
        x = x.T
    Sigmainv = np.identity(n)*(1/(h**2))  
    c,u = C(data, x, h) 
    const = constant(h, n, N) 


    H = (np.dot(c*u,u.T) - np.sum(c)*Sigmainv)/N
    return H
    

def C(data, x, h):
    "Calculate the u = (x-x_i)/h and c = k(||(x-x_i)/h||^2)"
    data = np.atleast_2d(data)
    n, N = data.shape
    x = np.atleast_2d(x)
    m, M = x.shape
    if M != 1:
        x = x.T
        x = np.reshape(x, (n, 1))
    u = (x - data)/(h)
    u2 = np.sum(u*u, axis=0)
    c = kern(u2)
    return c,u
    


def kern(x):
    "Gaussian Kernel Profile"
    k_x = np.exp(-x/2.0)
    return k_x 

    
    