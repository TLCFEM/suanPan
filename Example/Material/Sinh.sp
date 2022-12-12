node 1 0 0
node 2 1 0
node 3 2 0
node 4 3 0

material Tanh1D 1 1000
material Elastic1D 2 100
material Sinh1D 3 1000
material Viscosity01 4 1 5

element Spring01 1 1 2 1
element Spring01 2 2 3 2
element Spring01 3 3 4 3
element Damper01 4 2 3 4

mass 5 2 1 1
mass 6 3 1 1

fix2 1 1 1 4
fix2 2 2 1 2 3 4

amplitude Tabular 1 h

cload 1 1 1 1 3

initial velocity 100 1 2

hdf5recorder 1 Node U1 2 3
hdf5recorder 2 Node V1 2 3

step dynamic 1 1
set ini_step_size 5E-3
set fixed_step_size 1

converger RelIncreDisp 1 1E-11 10 1

analyze

# Node 2:
# Coordinate:
# 1.0000 0
# Displacement:
# 0.5855 0
# Resistance:
# 5.9921e+02 0
# Velocity:
# 22.0072 0
# Acceleration:
# -5.9921e+02 0
# 
# Node 3:
# Coordinate:
# 2.0000 0
# Displacement:
# 0.5389 0
# Resistance:
# 4.9280e+02 0
# Velocity:
# 8.4285 0
# Acceleration:
# -4.9180e+02 0
peek node 2 3

# save recorder 1 2

exit
