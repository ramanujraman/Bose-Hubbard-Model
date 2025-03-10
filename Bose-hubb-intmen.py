from qutip import *
import numpy as np
import matplotlib.pyplot as plt

# Parameters for Bose Hubbard Model
t= 1              # hopping strength
U= 10             # Potential interaction
mu=3              # Chemical Potential
N=3               # Number of Bosons
M=5               # Number of Lattice Sites

# Create the Hamiltonian in parts
# Chemical Potential Term
h_mu = sum([tensor([num(N) if i==j else qeye(N) for j in range(M)]) for i in range(M)])

# Interaction Term
site_num_ops = [tensor([num(N) if i==j else qeye(N) for j in range(M)]) for i in range(M)]
h_U = sum([site_num_ops[i]*(site_num_ops[i]-tensor([qeye(N) for j in range(M)])) for i in range(M)])

# Hopping Term
site_hop_ops = [tensor([create(N) if i==j else destroy(N) if i==j+1 else qeye(N) for j in range(M)]) for i in range(M)]
h_t =sum(site_hop_ops)
h_t = h_t + h_t.dag()

#Intermediate Regime 
def H(t,U,mu):
     return -t*h_t +(U/2)*h_U - mu*h_mu

# Variation of Energy with U
U_list = np.linspace(0,100,100)
E = []
for u in U_list:
    E.append(H(t,u,mu).groundstate()[0])
plt.plot(U_list,E)
plt.xlabel('U')
plt.ylabel('Energy')     # Energy vs U plot
plt.show()

# Variation of Energy with mu
mu_list = np.linspace(0,10,100)
E = []
for m in mu_list:
    E.append(H(t,U,m).groundstate()[0])
plt.plot(mu_list,E)
plt.xlabel('mu')
plt.ylabel('Energy')     # Energy vs mu plot
plt.show()

# Variation of Energy with t
t_list = np.linspace(0,10,100)
E = []
for i in t_list:
    E.append(H(i,U,mu).groundstate()[0])
plt.plot(t_list,E)
plt.xlabel('t')
plt.ylabel('Energy')     # Energy vs t plot
plt.show()