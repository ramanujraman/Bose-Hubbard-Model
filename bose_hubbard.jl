using QuantumToolbox
N=20
ω = 1.0
a = destroy(N)
H = ω*a'*a
println(eigenenergies(H))








