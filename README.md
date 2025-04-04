
# Certifying Asymmetry in the Configuration of Three Qubits

This repository contains our codes to reproduce the numerical and experimental results from the paper "Certifying Asymmetry in the Configuration of Three Qubits" by A. Taoutioui, G. Drótos, and T. Vértesi. It demonstrates how to:

1. Compute the mirror symmetry bound  $Q_{\text{mirror}}$
2. Identify the most asymmetric configuration using the biased $I_3$ witness
3. Perform experimental asymmetry certification on IBM Quantum devices

## Parametrization of Target Bloch Vectors

Let us choose the following parametrization for our three target Bloch vectors:

$$
\begin{align*}
&\vec{n}_1 = (0,0,1),\\
&\vec{n}_2 = (\sin{\theta_2},0,\cos{\theta_2}),\\
&\vec{n}_3 = (\sin{\theta_3}\cos{\phi_3},\sin{\theta_3}\sin{\phi_3},\cos{\theta_3}).
\end{align*}
$$

In general, our target trine states are defined by the set of angles { $\alpha_{12}$, $\alpha_{13}$, $\alpha_{23}$ } regardless of the chosen parametrization. These angles are determined by the dot products of the vectors $\vec{n}_1$, $\vec{n}_2 $, and $\vec{n}_3$:

$$\cos(\alpha_{12}) = \vec{n}_1 \cdot \vec{n}_2 = \cos{\theta_2},$$

$$\cos(\alpha_{13}) = \vec{n}_1 \cdot \vec{n}_3 = \cos{\theta_3},$$

$$\hspace{4.8cm}\cos(\alpha_{23}) = \vec{n}_2 \cdot \vec{n}_3 = \sin{\theta_2}\sin{\theta_3}\cos{\phi_3}+\cos{\theta_2}\cos{\theta_3}.$$


When $\phi_3 = k\pi$ where $k \in \mathbb{Z}$, our target Bloch vectors are coplanar; otherwise, they are non-coplanar.


## Mirror Symmetry Bound Calculation

The function `MirrorBound` from `utils.py` takes as arguments the angles $\alpha_{12}$, $\alpha_{13}$ and $\alpha_{23}$  and returns the bound $Q_{\text{mirror}}$. This function solves the optimization problem (29) subject to the mirror symmetry constraint in Eq. (30) along with the constraints on the preparation and measurement Bloch vectors (i.e., their norm can be less than or equals to 1). The problem is solved for the three different possible mirror symmetry constraints, yielding three bounds. The bound $Q_{\text{mirror}}$ is simply the maximum of these three bounds. 
The global optimization has been performed by using the Squential Least Squares Quadratic Programming  along with the basinhopping algorithm to avoid local maxima.

The mirror symmetry bound can be computed using the following piece of code:

```python
from utils import Degree, MirrorBound,  OverallBound

# Target angles: we take as example the most asymmetric configuration
alpha12 = 58.4 * Degree
alpha13 = 121.6 * Degree
alpha23 = 180 * Degree

Qmirror,[Qmirror123, Qmirror213, Qmirror312] = MirrorBound(alpha12, alpha13, alpha23)

print("Qmirror=", Qmirror)
# If you are intrested in the mirror bounds Q^{ijk}_{\text{mirror}} in Eq. ()
# uncomment the following lines
#print("Qmirror123=", Qmirror123)
#print("Qmirror213=", Qmirror213)
#print("Qmirror312=", Qmirror312)
````
which outputs `Qmirror= 5.82842712`. This obtained bound is tight, as it matches the same bound obtained when using the Lasserre hierarchy approach to relax our Quadratically Constrained Quadratic Programming (QCQP) problem into a Semidefinite Program (SDP). This task has been performed by using the Matlab script `Lassere_SDP_Q.m`.

Based on our numerical experience, the mirror bound can always be achieved using pure states. We have aslo included scripts in this repository that perform the optimization over pure states (see script names ending with `_pure`) and over pure coplanar states (see script names ending with `_coplanar_pure`). 

## Identifying the Most Asymmetric Configuration

The code below demonstrates how to use `utils.py` to identify the most asymmetric configuration by optimizing the gap in Eq. (32). The optimization has been performed by using the heuristic Nelder-Mead method. To ensure that the obtained maximum gap corresponds to a global maximum, we solve the optimization problem multiple times with random initial guesses (i.e. random simplices) and select the largest maximum as the global one. According to our tests, considering 10 random simplices is enough to find the global maximum. In this case, the execution time is approximately 25 minutes on a personal laptop.

```python
import numpy as np
from scipy.optimize import minimize
from utils import Degree, gap_witness

# Objective function (theta2, theta3, phi3):
def objective(x):
    alpha12, alpha13, phi3 = x  # alpha12 = theta2 and alpha13 = theta3
    # Calculation of alpha23
    alpha23 = np.arccos(np.sin(alpha12) * np.sin(alpha13) * np.cos(phi3)\
              + np.cos(alpha12) * np.cos(alpha13))
    return -gap_witness(alpha12, alpha13, alpha23)

# Heuristic optimization using Nelder-Mead method
num_initial_guesses = 10
best_cost = np.inf
best_result = None

for _ in range(num_initial_guesses):
    # Initial guess
    x0 = np.pi * np.random.rand(2)
    x0 = np.append(x0, 2 * np.pi * np.random.rand())

    # Solve using Nelder-Mead algorithm
    result = minimize(objective, x0, method='nelder-mead',\
             bounds=[(0.0, np.pi), (0.0, np.pi), (0, 2 * np.pi)], tol=1e-12)

    if result.fun < best_cost:
        best_cost = result.fun
        best_result = result

print("Delta_max = ",-best_result.fun)

alpha12 = best_result.x[0]
alpha13 = best_result.x[1]
phi3 = best_result.x[2]

alpha23 = np.arccos(np.sin(alpha12) * np.sin(alpha13) * np.cos(phi3) + np.cos(alpha12) * np.cos(alpha13))
print(alpha12 / Degree, alpha13 / Degree, alpha23 / Degree)
```
This script returns $\alpha_{12}=58.4^\circ$, $\alpha_{13}= 121.6^\circ$ and $\alpha_{23}= 180^\circ$ which define the most asymmetric configuration with the gap $\Delta_{max}\approx0.1111$.

## Experimental Asymmetry Certification
To experimentally certify the asymmetry property of a given target trine states, the script `IBMQuantumExperiment.py` can be used. One simply needs to specify the set of angles { $\alpha_{12}$, $\alpha_{13}$, $\alpha_{23}$ }, the number of shots to be considered, the IBM Quantum device (i.e., backend name), and the physical qubit number on which to run the designed quantum circuits.

After executing the circuits, the script `ExpDataProcessing.py` can be used to retrieve the jobs from the backend and evaluate the experimental value of the witness $I_6$, along with the corresponding standard deviation. If the obtained value violates the mirror bound $Q_{\text{mirror}}$ of the target trine states, the certification task is considered successful.

It is important to note that `IBMQuantumExperiment.py`, `ExpDataProcessing.py`  and `utils.py` scripts must be in the same directory, as `IBMQuantumExperiment.py` generates an input file for `ExpDataProcessing.py`. Before running the experiment, we select the best quality qubit based on the calibration data of the desired backend by specifying its number in the initial layout. This increases the likelihood of certifying the asymmetry property of our target. 

The `IBMQuantumExperiment.py` script is designed for coplanar targets (i.e., $\alpha_{12}+\alpha_{13}+\alpha_{23}=2\pi$). However, the script can be slightly modified to accommodate non-coplanar configurations using the aforementioned parametrization.

__Note__: Before using the above scripts, ensure that the necessary packages listed in `Requirements.txt` are installed. If not, run the following command to install them: `pip install -r Requirements.txt`.
