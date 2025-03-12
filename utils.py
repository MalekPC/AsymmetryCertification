from numpy import pi, sin, cos, sqrt
import numpy as  np
from scipy.optimize import minimize

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


def sample_a_Bloch_Vector():
    '''This function samples random a Bloch vector'''
    vec = np.random.normal(0, 1, 3)       # Sample from normal distribution
    vec /= np.linalg.norm(vec)            # Normalize to unit sphere
    r = np.random.uniform(0, 1) ** (1/2)  # sample the radius uniformly between 0 and 1  
    return vec * r


def I6(w12,w13,w23,x):
    ''' This is the witness function I6 with auxilary Bloch vectors chosen as: n4=- 
        (n1+n2)/norm(n1+n2), n5=-(n1+n3)/norm(n1+n3) and n6=-(n2+n3)/norm(n2+n3) (this is the 
        optimal choice of the auxilary Bloch vectors n4, n5, n6)
        w12, w13, w23 are the biases which are calculated from the relative angles between the 
        target pair of Bloch vectors.
    '''
    x1,y1,z1,\
    x2,y2,z2,\
    x3,y3,z3,\
    xm1,ym1,zm1,\
    xm2,ym2,zm2,\
    xm3,ym3,zm3,\
    xm4,ym4,zm4,\
    xm5,ym5,zm5,\
    xm6,ym6,zm6 = x # cordinates of the Bloch vectors of the preparations and measurements

    #preparations' Bloch vector  
    n1= np.array([x1,y1,z1])
    n2= np.array([x2,y2,z2])
    n3= np.array([x3,y3,z3])

    #measurements' Bloch vector 
    m1= np.array([xm1,ym1,zm1])
    m2= np.array([xm2,ym2,zm2])
    m3= np.array([xm3,ym3,zm3])
    m4= np.array([xm4,ym4,zm4])
    m5= np.array([xm5,ym5,zm5])
    m6= np.array([xm6,ym6,zm6])

    return  w12 + w13 + w23 \
        +w12*(n1+n2)@m1 + (1-w12)*(n1-n2)@m2 \
        +w13*(n1+n3)@m3 + (1-w13)*(n1-n3)@m4 \
        +w23*(n2+n3)@m5 + (1-w23)*(n2-n3)@m6


#========== Constraints:
ConstraintsOnBlochVectors=[
         {'type': 'ineq', 'fun': lambda x: 1 - (x[0]**2 + x[1]**2+x[2]**2)},            # |n1| ≤ 1 
         {'type': 'ineq', 'fun': lambda x: 1 - (x[3]**2 + x[4]**2+x[5]**2)},            # |n2| ≤ 1 
         {'type': 'ineq', 'fun': lambda x: 1 - (x[6]**2 + x[7]**2+x[8]**2)},            # |n3| ≤ 1 
         {'type': 'ineq', 'fun': lambda x: 1 - (x[9]**2 + x[10]**2+x[11]**2)},          # |m1| ≤ 1 
         {'type': 'ineq', 'fun': lambda x: 1 - (x[12]**2 + x[13]**2+x[14]**2)},         # |m2| ≤ 1 
         {'type': 'ineq', 'fun': lambda x: 1 - (x[15]**2 + x[16]**2+x[17]**2)},         # |m3| ≤ 1 
         {'type': 'ineq', 'fun': lambda x: 1 - (x[18]**2 + x[19]**2+x[20]**2)},         # |m4| ≤ 1 
         {'type': 'ineq', 'fun': lambda x: 1 - (x[21]**2 + x[22]**2+x[23]**2)},         # |m5| ≤ 1 
         {'type': 'ineq', 'fun': lambda x: 1 - (x[24]**2 + x[25]**2+x[26]**2)}          # |m6| ≤ 1 
]

# Mirror symmetry constraints [Eq. (7) in the manuscript:
#
# |n1-n2|=|n1-n3|:
def Mirror123Constraint(x): 
    return (x[0]-x[3])**2 + (x[1]-x[4])**2 + (x[2]-x[5])**2 -((x[0]-x[6])**2+(x[1]-x[7])**2+(x[2]-x[8])**2)
#
#|n2-n1|=|n2-n3|:
def Mirror213Constraint(x): 
    return (x[3]-x[0])**2 + (x[4]-x[1])**2 + (x[5]-x[2])**2 -((x[3]-x[6])**2+(x[4]-x[7])**2+(x[5]-x[8])**2)
#
#|n3-n1|=|n3-n2|:
def Mirror312Constraint(x): 
    return (x[6]-x[0])**2 + (x[7]-x[1])**2 + (x[8]-x[2])**2 -((x[6]-x[3])**2+(x[7]-x[4])**2+(x[8]-x[5])**2)
#================

def MirrorBound(alpha12, alpha13, alpha23):

    w12 = weight_function(alpha12)
    w13 = weight_function(alpha13)
    w23 = weight_function(alpha23)


    def objective(x):
        return -I6(w12,w13,w23,x)

    num_initial_guesses = 5
    best_cost = np.inf
    best_result = None

    for _ in range(num_initial_guesses):

         # Initial guess
         x0 = []
         for i in range(9):
             x0.append(sample_a_Bloch_Vector())

         initial_guess = np.concatenate(x0).tolist()

         # constraints         
         constraints1 = ConstraintsOnBlochVectors.copy() \
         + [{'type':   'eq', 'fun': Mirror123Constraint }]   # |n1-n2|=|n1-n3| 
         constraints2 = ConstraintsOnBlochVectors.copy() \
         + [{'type':   'eq', 'fun': Mirror213Constraint }]   # |n2-n1|=|n2-n3|
         constraints3 = ConstraintsOnBlochVectors.copy() \
         + [{'type':   'eq', 'fun': Mirror312Constraint }]   # |n3-n1|=|n3-n2|

         # Solve the optimization problem using Squential Least Squares Quadratic Programming  
         result1 = minimize(objective, initial_guess, method='SLSQP',\
         constraints=constraints1, tol=1e-10)
         result2 = minimize(objective, initial_guess, method='SLSQP',\
         constraints=constraints2, tol=1e-10)
         result3 = minimize(objective, initial_guess, method='SLSQP',\
         constraints=constraints3, tol=1e-10)

         result = min(result1.fun,result2.fun,result3.fun)

         if  result < best_cost:
             best_cost = result
             best_result = result
             #OptimalResults=[result1,result2,result3]

         Qmirror = -best_result
         return Qmirror    
         
# Gap witness in Eq. (32) in the manuscript:
def gap_witness(alpha12,alpha13,alpha23):
    return OverallBound(alpha12, alpha13, alpha23)-MirrorBound(alpha12, alpha13, alpha23)
