from numpy import pi, sin, cos, sqrt
import numpy as  np
from scipy.optimize import minimize, basinhopping

Degree= pi/180.0

def weight_function(alpha):
            return 0.5*(1.0 + (1.0 - abs(np.sin(alpha)))*(1.0/np.cos(alpha)))

def OverallBound(alpha12, alpha13, alpha23):
    w12 = weight_function(alpha12)
    w13 = weight_function(alpha13)
    w23 = weight_function(alpha23)
    """ See Eq. (27) in the manuscript """
    I12 = 2.0*np.sqrt(w12**2+(1.0-w12)**2) + w12
    I13 = 2.0*np.sqrt(w13**2+(1.0-w13)**2) + w13
    I23 = 2.0*np.sqrt(w23**2+(1.0-w23)**2) + w23
    return I12 + I13 + I23

def I6(z,w12,w13,w23,x):
    
    n = []
    for i in range(6):
    	n.append(np.array([0, 0, 0])) # initializing n
    n[z] = np.array([0, 0, 1])
    n[(z+1)%3] = np.array([sin(x[0])*cos(x[1]), sin(x[0])*sin(x[1]), cos(x[0])])
    n[(z+2)%3] = np.array([sin(x[0])*cos(x[1]), -sin(x[0])*sin(x[1]), cos(x[0])])
    n[3] = np.array([sin(x[2])*cos(x[3]), sin(x[2])*sin(x[3]), cos(x[2])])
    n[4] = np.array([sin(x[4])*cos(x[5]), sin(x[4])*sin(x[5]), cos(x[4])])
    n[5] = np.array([sin(x[6])*cos(x[7]), sin(x[6])*sin(x[7]), cos(x[6])])
    
    m = []
    m.append(w12*(n[0]+n[1]-n[3]))
    m.append((1-w12)*(n[0]-n[1]))
    m.append(w13*(n[0]+n[2]-n[4]))
    m.append((1-w13)*(n[0]-n[2]))
    m.append(w23*(n[1]+n[2]-n[5]))
    m.append((1-w23)*(n[1]-n[2]))
    for i in range(6):
	    m[i] /= np.linalg.norm(m[i])

    return  -w12*n[3]@m[0] - w13*n[4]@m[2] - w23*n[5]@m[4] \
        +w12*(n[0]+n[1])@m[0] + (1-w12)*(n[0]-n[1])@m[1] \
        +w13*(n[0]+n[2])@m[2] + (1-w13)*(n[0]-n[2])@m[3] \
        +w23*(n[1]+n[2])@m[4] + (1-w23)*(n[1]-n[2])@m[5]

def MirrorBound(alpha12, alpha13, alpha23):

    w12 = weight_function(alpha12)
    w13 = weight_function(alpha13)
    w23 = weight_function(alpha23)


    def objective1(x):
        return -I6(0,w12,w13,w23,x)
    def objective2(x):
        return -I6(1,w12,w13,w23,x)
    def objective3(x):
        return -I6(2,w12,w13,w23,x)

    # Initial guess
    initial_guess = []
    for i in range(4):
        initial_guess.append(np.random.uniform(0, pi))
        initial_guess.append(np.random.uniform(0, 2*pi))

    local_minimizer= {"method": "SLSQP"}
    
    result1= basinhopping(objective1, initial_guess, niter=100, minimizer_kwargs=local_minimizer)
    result2= basinhopping(objective2, initial_guess, niter=100, minimizer_kwargs=local_minimizer)
    result3= basinhopping(objective3, initial_guess, niter=100, minimizer_kwargs=local_minimizer)
   
    result = min(result1.fun,result2.fun,result3.fun)
    Qmirror = -result
    OptimalResults=[-result1.fun,-result2.fun,-result3.fun]
    return [Qmirror, OptimalResults]    
