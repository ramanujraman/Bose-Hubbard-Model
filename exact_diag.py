from qutip import *
import numpy as np
import matplotlib.pyplot as plt

# Fixed Parameters
U = 1              # On-site interaction strength
mu = 1             # Chemical potential
N = 5              # Local Hilbert space dimension (max bosons per site)
M = 4              # Number of lattice sites

# Create operators
num_ops = [tensor([num(N) if i == j else qeye(N) for j in range(M)]) for i in range(M)]
create_ops = [tensor([create(N) if i == j else qeye(N) for j in range(M)]) for i in range(M)]
destroy_ops = [tensor([destroy(N) if i == j else qeye(N) for j in range(M)]) for i in range(M)]

# Build static parts of the Hamiltonian (independent of t)
# On-site interaction: U/2 * n_i(n_i - 1)
H_U = sum([U / 2 * num_ops[i] * (num_ops[i] - 1) for i in range(M)])
# Chemical potential: -mu * n_i
H_mu = sum([-mu * num_ops[i] for i in range(M)])

# Function to construct the full Hamiltonian given t
def bose_hubbard_H(t):
    H_t = sum([create_ops[i] * destroy_ops[(i + 1) % M] + destroy_ops[i] * create_ops[(i + 1) % M] for i in range(M)])
    return -t * H_t + H_U + H_mu

# Vary hopping strength t and compute average boson number at site 2
t_list = np.linspace(0, 10, 50)
n_site2 = []

for t_val in t_list:
    H = bose_hubbard_H(t_val)
    gs = H.groundstate()[1]
    n2 = expect(num_ops[1], gs)  # Average boson number at site 2
    n_site2.append(n2)

# Plotting
plt.plot(t_list, n_site2, 'b-')
plt.xlabel('Hopping strength (t)')
plt.ylabel('Average boson number at site 2')
plt.title('n vs t in Bose-Hubbard Model')
plt.grid(True)
plt.show()
