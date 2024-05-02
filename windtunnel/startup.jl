using MPI

PROCS = [4]
FILE = "wind.jl"
timings = [;]

for (i, v) in enumerate(PROCS)
    println("procs: $v")
    t = @elapsed mpiexec(cmd->run(`$cmd -np $v julia --project=. $FILE`));
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