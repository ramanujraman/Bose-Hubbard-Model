from qutip import *
import numpy as np
import matplotlib.pyplot as plt

# Parameters for Bose Hubbard Model
t= 1            # hopping strength
U= 100             # Potential interaction
mu=3             # Chemical Potential
N=5             # Number of Bosons
M=5              # Number of Lattice Sites

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

print(H(t,U,mu).groundstate()[0])