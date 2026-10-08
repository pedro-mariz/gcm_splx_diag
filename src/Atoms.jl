# Superposition S (ss)
S_int(αi, αj) = ( π / (αi + αj) ) ^ 1.5

# Kinetic term T (ss)
T_int(αi, α2) = 3.0 * (αi * αj) / (αi + αj) * S_int(αi, αj)

# Nuclear attraction Vn (ss)
Vn_int(αi, αj; Z=2.0) = - 2.0 * π * Z / (αi + αj)

# Electron repulsion Vee (ssss)
Vee_int(αi, αj, αk, αl) = 2.0 * π ^2.5 / (αi + αj) / (αk + αl) / sqrt(αi + αj + αk + αl)

# Hamiltonian for single electron
h_int(αi, αj; Z=2.0) = T_int(αi, αj) + Vn_int(αi, αj; Z=Z)