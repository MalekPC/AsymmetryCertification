from utils_coplanar import Degree, MirrorBound,  OverallBound

# Target angles: we take as example the most asymmetric configuration
alpha12 = 54 * Degree
alpha13 = 112 * Degree
alpha23 = 194 * Degree

Qmirror,[Qmirror123,Qmirror213, Qmirror312] = MirrorBound(alpha12, alpha13, alpha23)
Qmax = OverallBound(alpha12, alpha13, alpha23)



print("Qmirror=", Qmirror)
print("Qmirror123=", Qmirror123)
print("Qmirror213=", Qmirror213)
print("Qmirror312=", Qmirror312)
print("Qmax=", Qmax)

print("delta=", Qmax-Qmirror)
