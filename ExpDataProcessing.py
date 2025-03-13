from qiskit_ibm_runtime import QiskitRuntimeService
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from numpy import pi, sin, cos, sqrt
import numpy as  np
from utils import weight_function, Degree, MirrorBound

# Reading the input parameters and the job ID from the external file
file_name = "ExperimentINPUT.txt"
with open(file_name, "r") as file:
    lines = file.readlines()
    alpha12 = float(lines[0].strip())
    alpha13 = float(lines[1].strip())
    alpha23 = float(lines[2].strip())
    shots = int(lines[3].strip())
    job_id = lines[4].strip()

print("alpha12:", round(alpha12/Degree, 8), "degrees")
print("alpha13:", round(alpha13/Degree, 8), "degrees")
print("alpha23:", round(alpha23/Degree, 8), "degrees")
print("Number of shots:", shots)
print("Job ID:", job_id)

# Choosing the service 

service = QiskitRuntimeService() 


def weight_function(alpha):
    ''' Calculating the bias parameter \( omega\) from the relative angle between
    the target pair of Bloch vectors \( \cos(\alpha) = \vec{n}_i \dot \vec{n}_j  \)'''
    return 0.5*(1.0 + (1.0 - abs(sin(alpha)))*(1./cos(alpha)))

#service = QiskitRuntimeService()
job = service.job(job_id)
#print("Job status:", job.status())

#  Retrieving the experiment counts:
result = job.result()
counts=[]
num_of_circuits = 3*2*3 # 3 qubits, 2 Pauli observables, 3 pairs of qubits

for  i in range(num_of_circuits):
        dist = result[i].data.c.get_counts()
        counts.append(dist)

# Calculating the associated expectation values
E=[]

#List of expectation values
for count in counts :
    evs=(count["0"]-count["1"])/shots
    E.append(evs)

# List (of dictionaries) of outcome probabilities 
proba=counts

for i in range(num_of_circuits):
    proba[i]["0"] = proba[i]["0"]/shots;
    proba[i]["1"] = proba[i]["1"]/shots

# Variances of the expectation values 
Var=[]
for i in range(num_of_circuits):
    Var.append(4*proba[i]["0"]*(1-proba[i]["0"])/shots)

# Calculating witness I3(\omgea_{12})

Exy =  np.array([[E[0], E[1]],
                 [E[2], E[3]],
                 [E[4], E[5]]] );

Varxy = np.array([[Var[0],  Var[1]],
                  [Var[2],  Var[3]],
                  [Var[4],  Var[5]]] );

w=weight_function(alpha12)

Wxy= np.array([[w,   (1-w)],
               [w,  -(1-w)],
               [-w,  0]]);

I_12=np.sum(Exy*Wxy)

# standard deviation 
Var12 = np.sum(Wxy**2*Varxy)
delta12=sqrt(Var12)

print("Witness I_3(omega12) = ", round(I_12,6), "±", round(delta12,6))

# Calculating witness I3(\omgea_{13})

Exy =  np.array([[E[6], E[7]],
                 [E[8], E[9]],
                 [E[10], E[11]]] );

Varxy = np.array([[Var[6],  Var[7]],
                  [Var[8],  Var[9]],
                  [Var[10], Var[11]]] );

w=weight_function(alpha13)

Wxy= np.array([[w,   (1-w)],
               [w,  -(1-w)],
               [-w,  0]]);

I_13=np.sum(Exy*Wxy)


# standard deviation 
Var13 = np.sum(Wxy**2*Varxy)
delta13=sqrt(Var13)

print("Witness I_3(omega13) = ", round(I_13,6), "±", round(delta13,6))

# Calculating witness I3(\omgea_{23})

Exy =  np.array([[E[12], E[13]],
                 [E[14], E[15]],
                 [E[16], E[17]]] );

Varxy = np.array([[Var[12],  Var[13]],
                  [Var[14],  Var[15]],
                  [Var[16],  Var[17]]] );

w=weight_function(alpha23)

Wxy= np.array([[w,   (1-w)],
               [w,  -(1-w)],
               [-w,  0]]);

I_23=np.sum(Exy*Wxy)


# standard deviation 
Var23 = np.sum(Wxy**2*Varxy)
delta23=sqrt(Var23)

print("Witness I_3(omega23) = ", round(I_23,6), "±", round(delta23,6))

# Overall witness
I_6 = I_12 + I_13 + I_23

# standard deviation
delta = sqrt( delta12**2+delta13**2+delta23**2)

print("Witness I_6 = ", round(I_6,6), "±" , round(delta,6))

print("Mirror symmetry bound: ", round(MirrorBound(alpha12, alpha13, alpha23),6))
