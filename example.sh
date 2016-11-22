
# update the path to fastpm

mpirun -n 4 ~/source/fastpm/src/fastpm bluetides-pathfinder-50.lua

mpirun -n 4 python preion-make-zreion.py bluetides-pathfinder-50/fastpm_0.0909 bluetides-pathfinder-50/UVModulation_0.0909
