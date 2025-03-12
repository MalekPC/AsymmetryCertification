from qiskit import QuantumRegister, QuantumCircuit
import numpy as np
from numpy import pi,cos, sin, sqrt, exp, arccos
from scipy.linalg import norm
from qiskit_ibm_runtime import QiskitRuntimeService
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
# Running the transpiled quantum circuits 
from  qiskit_ibm_runtime import SamplerV2  as Sampler 
from  qiskit_ibm_runtime import SamplerOptions 


# Choosing the service 
token="Place your IBM's API token here before running the code"
service = QiskitRuntimeService(channel="ibm_quantum", token=token) 


# Angles defining the asymmetric target trine states: 
alpha12 = 58.4*Degree    #  cos(alpha12) = \vec{n}_1 \dot  \vec{n}_2
alpha13 = 121.6*Degree   #  cos(alpha13) = \vec{n}_1 \dot  \vec{n}_3
alpha23 = 180*Degree     #  cos(alpha23) = \vec{n}_2 \dot  \vec{n}_3

# Number of shots
shots=4096

# Choosing the IBM Quantum backend
backend_name = "ibm_brisbane"

#Qubit number (before running the code we choose the best quality qubit based on the backend calobraion data provided by IBM Quantum Cloud).

qubit_number=9
#============================ Functions ===================================== 
''' Here we define some useful fucntions which will be used later on to build 
the quantum circuits intented to certify our target triple of qubit states.'''

def cartesian_to_spherical(x,y,z):
    XsqPlusYsq = x**2 + y**2
    r = sqrt(XsqPlusYsq + z**2)                      # r
    theta= np.arctan2(sqrt(XsqPlusYsq),z)            # theta
    phi =  np.arctan2(y,x)                           # phi
    return [theta, phi,r]

def weight_function(alpha):
    ''' Calculating the bias parameter \( omega\) from the relative angle between
    the target pair of Bloch vectors \( \cos(\alpha) = \vec{n}_i \dot \vec{n}_j  \)'''
    return 0.5*(1.0 + (1.0 - abs(sin(alpha)))*(1./cos(alpha)))
    
def projective_meas_gate(m,meas_label):
    '''This function is useful for measuring a Pauli observable represented by the unit Bloch vector \( \vec{m} \).
       It returns a single qubit gate that must be applied just before the measurement in the computational basis. 
       The resulting statistics will be equivalent to performing a measurement in the basis defined by the Pauli 
       matrix associated with the Bloch vector \( \vec{u} \).'''
    x,y,z=m
    theta, phi,r = cartesian_to_spherical(x,y,z)
    
    if round(r,10) !=1:
       raise ValueError("Error: the measurment Bloch vector must be unit vector!")

    name="m"+str(meas_label)
    
    """ Building the gate (i.e., unitary operator )"""
    u = np.array([[cos(theta/2.0), exp(-1j*phi)*sin(theta/2.0)],
                 [cos((pi-theta)/2.0), exp(-1j*(phi+pi))*sin((pi-theta)/2.0)]])    
    
    """------- Creating the measurement gate using the "unitary" function in Qiskit---"""
    meas_gate = QuantumCircuit(1,1) 
    meas_gate.unitary( u ,0, label=name)
    meas_gate.measure([0],[0])           
    return  meas_gate   

def MUB_pair(vector1, vector2):
    ''' This function returns two measurement Bloch vectors associated with the mutually unbiased baises (MUB) 
        defned by  the preparation Bloch vectors vector1 and vector2 '''  
    alpha=arccos(vector1@vector2)
    w=weight_function(alpha)
    wrounded=round(w,10)
    
    if  wrounded !=0.0 or  wrounded!= 1.0 : 
          u=sqrt(w**2+(1-w)**2)/(2.0*w)*(vector1+vector2)
          v=sqrt(w**2+(1-w)**2)/(2.0*(1-w))*(vector1-vector2)
    
    if  wrounded==0.0 : # --> alpha= pi
       v=vector1
       x0,y0,z0 = vector1
       Theta,Phi,R =cartesian_to_spherical(x0,y0,z0) 
       u=[sin(Theta+pi/2.0)*cos(Phi),sin(Theta+pi/2.0)*sin(Phi),cos(Theta+pi/2.0)]
        
    if  wrounded==1.0 : # --> alpha= 0
       u=vector1
       x0,y0,z0 = vector1
       Theta,Phi,R =cartesian_to_spherical(x0,y0,z0) 
       v=[sin(Theta+pi/2.0)*cos(Phi),sin(Theta+pi/2.0)*sin(Phi),cos(Theta+pi/2.0)]
    return [u,v]

def RY(vector, theta):
       """This function rotate in 3D vector around the Y-axis, 
          one need to specify the vector to be rotated and 
          the angle of rotation \(\theta\)"""

       U=np.array([[cos(theta),0.0,sin(theta)],
             [0,        1.0,        0.0],
             [-sin(theta),0.0,cos(theta)]
       ])
       return  U@vector

def build_circuit(prepared_Bloch_vector, measurement_Bloch_vector,meas_label):
    n=prepared_Bloch_vector
    m=measurement_Bloch_vector
    
    # preparation
    circuit=QuantumCircuit(1,1)
    x,y,z=n
    theta, phi,r= cartesian_to_spherical(x,y,z) 

    circuit.u(theta,phi,0, qubit=0)

    circuit.barrier() 

    #measurement
   
    meas=projective_meas_gate(m,meas_label)
    circuit=circuit.compose(meas)

    return circuit
#==========================================================================

#== Preparation Bloch vectors (target + auxilary Bloch vectors)

# Target
n1=np.array([0,0,1])
n2=RY(n1,alpha12)
n3=RY(n2,alpha23)

# Auxilary 
n4=-(n1+n2)/norm(n1+n2)
n5=-(n1+n3)/norm(n1+n3)
n6=-(n2+n3)/norm(n2+n3)

#== Building the optimal measurement Bloch vectors

m1,m2 = MUB_pair(n1,n2)

m3,m4 = MUB_pair(n1,n3)

m5,m6 = MUB_pair(n2,n3)

#=== Theoretical overall witness Q 

# biases
w12 = weight_function(alpha12)
w13 = weight_function(alpha13)
w23 = weight_function(alpha23)

I12=w12*(n1+n2-n4)@m1 +(1.0-w12)*(n1-n2)@m2
print("I_12=",round(I12,6))  

I13=w13*(n1+n3-n5)@m3 +(1.0-w13)*(n1-n3)@m4
print("I_13=",round(I13,6))  

I23= w23*(n2+n3-n6)@m5 +(1.0-w23)*(n2-n3)@m6
print("I_23=",round(I23,6))     

Q=(I12+I13+I23)
print("Q=",round(Q,6))

# Building the Quantum Circuits 

# Preparation and Measurement Bloch vectors
preparation1=[n1,n2,n4]
preparation2=[n1,n3,n5]
preparation3=[n2,n3,n6]
pair1=[m1,m2]
pair2=[m3,m4]
pair3=[m5,m6]

circuits=[]
circuits1=[] 
circuits2=[] 
circuits3=[]

for i in range(3):
      for j in range(2):
         #=============================
         # Building the quantum circuit
         circuit1=build_circuit(preparation1[i], pair1[j], meas_label=j+1)
         circuit2=build_circuit(preparation2[i], pair2[j], meas_label=j+1)
         circuit3=build_circuit(preparation3[i], pair3[j], meas_label=j+1)

         circuits1.append(circuit1)
         circuits2.append(circuit2)
         circuits3.append(circuit3)

'''## Notes
circuits1 -----> (n1,m1); (n1,m2); (n2,m1); (n2,m2); .....
...
''' 

circuits = circuits1 + circuits2 + circuits3


backend = service.backend(backend_name)

# Transpiling the quantum circuits
pass_manager = generate_preset_pass_manager(optimization_level=3, backend=backend, initial_layout=[qubit_number])
circuits_transpiled = pass_manager.run(circuits)

sampler = Sampler(backend)
sampler.options.default_shots = shots


## run the circuits
job = sampler.run(circuits_transpiled)

job_id = job.job_id()
print("Job ID :", job_id)

# Saving the input parameters and the job ID to an external file
file_name = "ExperimentINPUT.txt"
with open(file_name, "w") as file:
    file.write(str(alpha12) + "\n")
    file.write(str(alpha13) + "\n")
    file.write(str(alpha23) + "\n")
    file.write(str(shots) + "\n")
    file.write(job_id) 