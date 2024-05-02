using MPI

PRECISION = 10000
STEPS = 8

# Structure to represent a solid particle in the tunnel surface.
struct Particle
    position::Vector{Int} # As [column, row]
    mass::Int
    # resistance::Int
    resistance::Float16
    speed::Vector{Float16} # Two dimensional velocity vector.
    oldFlow::Int
end

# TODO: Function to get wall time.
# TODO: Function to simplify matrix/2D array access.
# TODO: Update flow in a matrix position.
# TODO: Move particle function.
# TODO: Print the current state of the simulation.
# TODO: Print usage line in stderr.
# TODO: Start global timer
# TODO: Parallel initialization.
# TODO: Simulation.
# TODO: Particle movement each STEPS iterations.
# TODO: Effects due to particles each STEPS iterations.
# TODO: Propagation stage.
# TODO: MPI: Fill result arrays used for later output.
# TODO: Stop global timer.
# TODO: Output for the leaderboard.

function printArguments(rows, cols, max_iter, threshold, fan_pos, fan_size,
    fixed_band_pos, fixed_band_size, fixed_density, 
    moving_band_pos, moving_band_size, moving_density, 
    random_seq, n_particles, particles)
    println("### Print Arguments ###\nRows: $rows, Columns: $cols, Max iter: $max_iter, Threshold: $threshold\nFan: $fan_pos, $fan_size\nFixed particles : $fixed_band_pos $fixed_band_size $fixed_density\nMoving particles: $moving_band_pos $moving_band_size $moving_density")
    for i in 1:3
        println("Random sequence: $(random_seq[i])")
    end
    print("Number of particles: $n_particles\nparticles: ")
    for i in 1:n_particles
        print("$(particles[i]) ")
    end
end

"""Given a size the function returns a the dimensions for the cartesian grid as [rows, columns]"""
function getDims(n_ranks)
    if n_ranks == 1
        return [1, 1]
    elseif n_ranks == 2
        return [1, 2]
    elseif n_ranks == 4
        return [2, 2]
    elseif n_ranks == 8
        return [4, 2]
    elseif n_ranks == 16
        return [4, 4]
    end
end

"""For DAS the arguments are passed as parameters to main().
For local execution the arguments are in the command line.

Order of arguments: rows, cols, max_iter, threshold, fan_pos, fan_size, 
fixed_band_pos, fixed_band_size, fixed_density, 
moving_band_pos, moving_band_size, moving_density, 
rand_seq_1, rand_seq_2, rand_seq_3, 
particle_row, particle_col, particle_density, ...
"""
function main()
    
    # Simulation data.
    rows = parse(Int, ARGS[1])          # Tunnel rows.
    cols = parse(Int, ARGS[2])          # Tunnel columns.
    max_iter = parse(Int, ARGS[3])      # Max number of sim (simulation) steps.
    threshold = parse(Float16, ARGS[4]) # Threshold of variability to continue the sim.
    fan_pos = parse(Int, ARGS[5])       # First position of the fan.
    fan_size = parse(Int, ARGS[6])      # Fan size.
    fixed_band_pos = parse(Int, ARGS[7])        # First position of the band where FIXED particles Start.
    fixed_band_size = parse(Int, ARGS[8])       # Size of the band where FIXED particles start.
    fixed_density = parse(Float16, ARGS[9])     # Density of starting FIXED particles.
    moving_band_pos = parse(Int, ARGS[10])      # First position of the band where MOVING particles start.
    moving_band_size = parse(Int, ARGS[11])     # Size of the band where MOVING particles start.
    moving_density = parse(Float16, ARGS[12])   # Density of starting MOVING particles.
    flow = []           # Tunnel air flow.
    flow_copy = []      # Tunnel air flow (ancillary, of secondary importance, copy).
    particle_locs = []  # Quickly locate particle positions.

    # Random sequences initializer.
    random_seq = [ARGS[i] for i in 13:15]

    # Initialize MPI basics.
    MPI.Init()
    comm = MPI.Comm_dup(MPI.COMM_WORLD)
    n_ranks = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)

    # Initialize MPI cartesian grid.
    dims = getDims(n_ranks)
    periodic = map(_ -> false, dims)
    reorder = false
    comm_cart = MPI.Cart_create(comm, dims; periodic, reorder)

    # Read particles from the arguments.
    n_particles = (length(ARGS) - 15) ÷ 3               # Number of particles.
    particles = Vector{Particle}(undef, n_particles)    # List to store cells information
    for i in 1:n_particles
        offset = i
        row = parse(Int, ARGS[13 + offset * 3]) * PRECISION
        col = parse(Int, ARGS[14 + offset * 3]) * PRECISION
        res = parse(Float16, ARGS[15 + offset * 3]) * PRECISION
        particles[i] = Particle([row, col], 0, res, [0,0], 0)
    end

    # TODO: Generate fixed/moving particles in the band.
    

    if rank == 0
        printArguments(rows, cols, max_iter, threshold, fan_pos, fan_size,
            fixed_band_pos, fixed_band_size, fixed_density, 
            moving_band_pos, moving_band_size, moving_density, random_seq, n_particles, particles)
    end
    MPI.Barrier(comm_cart)

    #=
    start global timer

    init.

    sim.
        1. change inlet values
        2. move particles
        3. effects due to particles
        4. copy into ancillary structure
        5. propagation stage
    =#

end

main()
