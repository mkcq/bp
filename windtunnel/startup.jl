using MPI

PROCS = [4]
FILE = "wind.jl"
timings = [;]

#=
rows, cols, max_iter, threshold, fan_pos, fan_size, 
fixed_band_pos, fixed_band_size, fixed_density, 
moving_band_pos, moving_band_size, moving_density, 
rand_seq_1, rand_seq_2, rand_seq_3, 
particle_row, particle_col, particle_density,
=#

rows = 32               # Tunnel rows.
cols = 32               # Tunnel columns.
max_iter = 1            # Max number of sim (simulation) steps.
threshold = 0.          # Threshold of variability to continue the sim.
fan_pos = 2             # First position of the fan.
fan_size = 28           # Fan size.
fixed_band_pos = 1      # First position of the band where FIXED particles Start.
fixed_band_size = 0     # Size of the band where FIXED particles start.
fixed_density = 0       # Density of starting FIXED particles.
moving_band_pos = 5     # First position of the band where MOVING particles start.
moving_band_size = 3    # Size of the band where MOVING particles start.
moving_density = 0.3
random_seq_1 = 3434
random_seq_2 = 1242
random_seq_3 = 965


for (i, v) in enumerate(PROCS)
    println("procs: $v")
    t = @elapsed mpiexec(cmd->run(`$cmd -np $v julia --project=. $FILE $rows $cols $max_iter $threshold $fan_pos $fan_size $fixed_band_pos $fixed_band_size $fixed_density $moving_band_pos $moving_band_size $moving_density $random_seq_1 $random_seq_2 $random_seq_3 8 17 0.9`));
    push!(timings, (v, t))
end
display(timings)

#=
#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=16

. /etc/profile.d/lmod.sh
module load openmpi/gcc/64
. /etc/profile.d/lmod.sh
module load julia/1.7.3

pwd

# Parallel

mpiexec -np 4 julia --project=. -e '
    include("wind.jl")
    rows = 1
    cols = 1
    max_iter = 1
    threshold = 1
    main(rows, cols, max_iter, threshold)
'
=#

#=
main(rows, cols, max_iter, var_threshold, inlet_pos, inlet_size, 
    particles_f_band_pos, particles_f_band_size, particles_f_density, 
    particles_m_band_pos, particles_m_band_size, particles_m_density,
    rand_seq_1, rand_seq_2, rand_seq_3, 
    particle_row, particle_col, particle_density, ...)
=#