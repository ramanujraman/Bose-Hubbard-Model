from qutip import *
import numpy as np
import matplotlib.pyplot as plt

# Parameters for Bose Hubbard Model
t=0.1            # hopping strength
U=2              # Potential interaction
mu=3             # Chemical Potential
N=3              # Number of Bosons
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

# Strong Interaction Limit (t=0)
# Total Hamiltonian for strongly interacting limit
H_strongly_interacting_limit = (U/2)*h_U - mu*h_mu
H_strongly_interacting_limit.groundstate()


# Weak Interaction Limit (U=0)
# Total Hamiltonian for weak interacting limit
H_weak_interacting_limit = -t*h_t - mu*h_mu
H_weak_interacting_limit.groundstate()

#Intermediate Regime 
H = -t*h_t +(U/2)*h_U - mu*h_mu
H.groundstate()