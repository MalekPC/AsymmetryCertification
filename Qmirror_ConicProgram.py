import cvxpy as cp
import numpy as np
from numpy import * 
from  scipy.optimize import minimize, basinhopping

Degree=np.pi/180

def weight_function(alpha):
    ''' given an angle alpha between two unit Bloch vectors, this function returns
    the associated weight omega (see Eq. (23) in manuscript)'''
    return 0.5*(1.0 + (1.0 - abs(np.sin(alpha)))*(1.0/np.cos(alpha)))


def mirror_bound(ijk, omega):
    ''' this function returns the bound Q^{ijk}_{max} by solving the conic convex program (B8) in the appendix B'''

    w12, w13, w23 = omega

    a1, b1 = w12, 1-w12
    a2, b2 = w13, 1-w13
    a3, b3 = w23, 1-w23

    # Toggle pure vs mixed
    pure = True #-> pure states (G_ii = 1). False -> mixed (0 <= G_ii <= 1)

    # Optional symmetry constraints (toggle as needed)
    sym123 = False
    sym213 = False
    sym312 = False

    if ijk == "123": sym123 = True
    if ijk == "213": sym213 = True
    if ijk == "312": sym312 = True

    # ------ Variables ------
    G = cp.Variable((3,3), PSD=True)     # Gram matrix (symmetric PSD)
    u12p = cp.Variable(nonneg=True)      # ||\vec{n}_1 + \vec{n}_2||
    u12m = cp.Variable(nonneg=True)      # ||\vec{n}_1 - \vec{n}_2||
    u13p = cp.Variable(nonneg=True)      # ||\vec{n}_1 + \vec{n}_3||
    u13m = cp.Variable(nonneg=True)      # ||\vec{n}_1 - \vec{n}_3||
    u23p = cp.Variable(nonneg=True)      # ||\vec{n}_2 + \vec{n}_3||
    u23m = cp.Variable(nonneg=True)      # ||\vec{n}_2 - \vec{n}_3||

    # Convenience indices
    r1, r2, r3    = G[0,0], G[1,1], G[2,2] 
    c12, c13, c23 = G[0,1], G[0,2], G[1,2]

    constraints = []

    # Diagonal bounds (mixed/pure)
    if pure:
        constraints += [r1 == 1, r2 == 1, r3 == 1]
    else:
        constraints += [r1 >= 0, r2 >= 0, r3 >= 0, 
                        r1 <= 1, r2 <= 1, r3 <= 1 ]
 
    # Symmetry options
    if sym123:
        constraints += [c12 == c13] 
    if sym213:
        constraints += [c12 == c23] 
    if sym312:
        constraints += [c13 == c23]

    # Square-root SOC constraints (accuracy-friendly)
    constraints += [
        u12p <= cp.sqrt(r1 + r2 + 2*c12),
        u12m <= cp.sqrt(r1 + r2 - 2*c12),
        u13p <= cp.sqrt(r1 + r3 + 2*c13),
        u13m <= cp.sqrt(r1 + r3 - 2*c13),
        u23p <= cp.sqrt(r2 + r3 + 2*c23),
        u23m <= cp.sqrt(r2 + r3 - 2*c23)
    ]

    #  Note: G is PSD by construction (Variable(..., PSD=True))

    # Objective: maximize linear form in u's
    objective0 = cp.Maximize(
        a1*u12p + b1*u12m +
        a2*u13p + b2*u13m +
        a3*u23p + b3*u23m)

    prob = cp.Problem(objective0, constraints)

    # Choose a solver; MOSEK is best if available.
    try:
       prob.solve(solver=cp.MOSEK,
               mosek_params={
                   "MSK_DPAR_INTPNT_TOL_REL_GAP": 1e-10,
                   "MSK_DPAR_INTPNT_CO_TOL_REL_GAP": 1e-10,
               },
               verbose=False)
    except Exception:
        prob.solve(solver=cp.MOSEK, eps=1e-12, max_iters=50000, verbose=False)


    return  prob.value


# Let us consider the following target:

theta2 = 60*Degree

n1 = np.array([0,0,1])  
n2 = np.array([sin(theta2),0,cos(theta2)])  
n3 = np.array([0,0,-1])  

# Compute the weights   
w12 = weight_function(np.arccos(np.dot(n1,n2)))
w13 = weight_function(np.arccos(np.dot(n1,n3)))
w23 = weight_function(np.arccos(np.dot(n2,n3)))
x = [w12, w13, w23]

# Compute the bounds  Q_^{ijk}_{max} 
Q123 = mirror_bound(ijk="123", omega=x)
Q213 = mirror_bound(ijk="213", omega=x)
Q312 = mirror_bound(ijk="312", omega=x)

# Take the max to evaluate the Q_mirror bound associated with the chosen target
Qmirror = np.max([Q123,Q213,Q312])

# Compute the overall witness (ie, witness imposing any of the symmetry constraints)
Qmax =  2*(sqrt(x[0]**2+(1-x[0])**2)+sqrt(x[1]**2+(1-x[1])**2)+sqrt(x[2]**2+(1-x[2])**2))

# Evaluate the gap between the Qmax and the Q_mirror
Gap = Qmax-Qmirror

print("Qmax=", np.round(Qmax,8))
print("Q_mirror=", np.round(Qmirror,8))
print("Associated Gap =", np.round(Gap, 8))
