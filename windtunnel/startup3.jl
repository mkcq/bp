using MPI

N_MACHINES = [4]
FILE = "wind3.jl"
TIMINGS = []

rows = 40
cols = 44
max_iter = 3710
threshold = 0.8
fan_pos = 2
fan_size = 40

# First position of the band where FIXED particles Start.
fixed_band_pos = 1
# Size of the band where FIXED particles start.
fixed_band_size = 0
# Density of starting FIXED particles.
fixed_density = 0

# First position of the band where MOVING particles start.
moving_band_pos = 5
# Size of the band where MOVING particles start.
moving_band_size = 3
# Density of starting MOVING particles.
moving_density = 0.3

random_seq_1 = 3434
random_seq_2 = 1242
random_seq_3 = 965

for (i, v) in enumerate(N_MACHINES)

    t = @elapsed mpiexec(cmd->run(`$cmd -np $v julia --project=. $FILE $rows $cols $max_iter $threshold $fan_pos $fan_size $fixed_band_pos $fixed_band_size $fixed_density $moving_band_pos $moving_band_size $moving_density $random_seq_1 $random_seq_2 $random_seq_3 8 17 0.9 9 16 0.9 9 17 0.9 9 18 0.9 10 16 0.9 10 17 0.9 10 18 0.9 11 16 0.9 11 17 0.9 11 18 0.9 12 16 0.9 12 17 0.9 12 18 0.9 13 16 0.9 13 17 0.9 13 18 0.9 14 7 0.9 14 8 0.9 14 9 0.9 14 10 0.9 14 11 0.9 14 12 0.9 14 13 0.9 14 14 0.9 14 15 0.9 14 16 0.9 14 17 0.9 14 18 0.9 14 19 0.9 14 20 0.9 14 21 0.9 14 22 0.9 14 23 0.9 14 24 0.9 14 25 0.9 14 26 0.9 14 27 0.9 15 6 0.9 15 7 0.9 15 8 0.9 15 9 0.9 15 10 0.9 15 11 0.9 15 12 0.9 15 13 0.9 15 14 0.9 15 15 0.9 15 16 0.9 15 17 0.9 15 18 0.9 15 19 0.9 15 20 0.9 15 21 0.9 15 22 0.9 15 23 0.9 15 24 0.9 15 25 0.9 15 26 0.9 15 27 0.9 15 28 0.9 16 5 0.9 16 6 0.9 16 7 0.9 16 8 0.9 16 9 0.9 16 10 0.9 16 11 0.9 16 12 0.9 16 13 0.9 16 14 0.9 16 15 0.9 16 16 0.9 16 17 0.9 16 18 0.9 16 19 0.9 16 20 0.9 16 21 0.9 16 22 0.9 16 23 0.9 16 24 0.9 16 25 0.9 16 26 0.9 16 27 0.9 16 28 0.9 16 29 0.9 17 5 0.9 17 6 0.9 17 7 0.9 17 8 0.9 17 9 0.9 17 10 0.9 17 11 0.9 17 12 0.9 17 13 0.9 17 14 0.9 17 15 0.9 17 16 0.9 17 17 0.9 17 18 0.9 17 19 0.9 17 20 0.9 17 21 0.9 17 22 0.9 17 23 0.9 17 24 0.9 17 25 0.9 17 26 0.9 17 27 0.9 17 28 0.9 17 29 0.9 18 6 0.9 18 7 0.9 18 8 0.9 18 9 0.9 18 10 0.9 18 11 0.9 18 12 0.9 18 13 0.9 18 14 0.9 18 15 0.9 18 16 0.9 18 17 0.9 18 18 0.9 18 19 0.9 18 20 0.9 18 21 0.9 18 22 0.9 18 23 0.9 18 24 0.9 18 25 0.9 18 26 0.9 18 27 0.9 18 28 0.9 19 16 0.9 19 17 0.9 19 18 0.9 20 16 0.9 20 17 0.9 20 18 0.9 21 16 0.9 21 17 0.9 21 18 0.9 22 16 0.9 22 17 0.9 22 18 0.9 23 16 0.9 23 17 0.9 23 18 0.9 24 16 0.9 24 17 0.9 24 18 0.9 25 16 0.9 25 17 0.9 25 18 0.9 26 16 0.9 26 17 0.9 26 18 0.9 27 16 0.9 27 17 0.9 27 18 0.9 28 17 0.9 29 17 0.9 30 17 0.9 31 16 0.9 31 17 0.9 31 18 0.9 32 15 0.9 32 16 0.9 32 17 0.9 32 18 0.9 32 19 0.9 33 15 0.9 33 16 0.9 33 17 0.9 33 18 0.9 33 19 0.9 34 17 0.9`));
    push!(TIMINGS, (v, t))

end

println("\n -> startup3.jl\nTimings = ")
display(TIMINGS)