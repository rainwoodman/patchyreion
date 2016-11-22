-- parameter file
------ Size of the simulation -------- 

-- For Testing
nc = 128
boxsize = 50.0

-------- Time Sequence ----
-- linspace: Uniform time steps in a
-- time_step = linspace(0.025, 1.0, 39)
-- logspace: Uniform time steps in loga
-- time_step = linspace(0.01, 1.0, 10)
time_step = linspace(0.0625, 1. / (10 + 1), 5)

output_redshifts= {10.0}  -- redshifts of output

-- Cosmology --
omega_m = 0.2814
h       = 0.697
sigma8  = 0.820
--
-- Start with a power spectrum file
-- Initial power spectrum: k P(k) in Mpc/h units
-- Must be compatible with the Cosmology parameter
read_powerspectrum= "powerspectrum-wmap9-linear.txt"
random_seed= 181170


pm_nc_factor = 2            -- Particle Mesh grid pm_nc_factor*nc per dimension in the beginning

np_alloc_factor= 4.0      -- Amount of memory allocated for particle

-------- Output ---------------

-- Dark matter particle outputs (all particles)
write_snapshot= "bluetides-pathfinder-50/fastpm" 
-- 1d power spectrum (raw), without shotnoise correction
write_powerspectrum = "bluetides-pathfinder-50/powerspec"

