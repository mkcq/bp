using MPI
using Random

PRECISION = 10000
STEPS = 8

# Structure to represent a solid particle in the tunnel surface.
mutable struct Particle
    position::Vector # As [column, row]
    mass::Int
    resistance::Int
    speed::Vector # Two dimensional velocity vector.
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
    println("Number of particles: $n_particles\nparticles:")
    for (k,v) in enumerate(particles)
        println("$k : $v")
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

function printDEBUG(n_fixed_band, n_moving_band, n_particles)
    x=length(ARGS)
    y=x-15
    z=y÷3
    println("$x, $y, $z\nn_fixed_band=$n_fixed_band\nn_moving_band=$n_moving_band\nn_particles=$n_particles")
end

function updateFan(iter, fan_pos, fan_size, flow)
    if iter % STEPS == 1
        for j in fan_pos:fan_pos+fan_size
            phase = iter ÷ STEPS * (π ÷ 4)
            phase_step = π ÷ 2 ÷ fan_size
            pressure_level = 9 + 2 * sin(phase + (j - fan_pos) * phase_step)
            noise = 0.5 - rand(1)[1]
            # Store level in the first row of the ancillary structure.
            flow[1, j] = trunc(Int, PRECISION * (pressure_level + noise))
        end
    end
end

function particleMovements(iter, rows, cols, particle_locs, particles, flow)
    if iter % STEPS == 1
        # Clean particle positions.
        for r in 1:iter
            if iter > rows
                break
            end
            for c in 1:cols
                particle_locs[r,c] = 0
            end
        end

        for p in 1:n_particles
            mass = particles[p].mass
            if mass == 0
                continue
            end
            moveParticle(flow, particles, p, rows, cols) # TODO: Implement this
        end

        # Annotate position.
        for p in 1:n_particles
            row = particles[p].position[1] ÷ PRECISION
            col = particles[p].position[2] ÷ PRECISION + 1
            particle_locs[row, col] += 1
        end
    end
end

function particleEffects(iter, cols, n_particles, particles, flow, flow_copy, particle_locs)
    if iter % STEPS == 1
        for p in 1:n_particles
            row = particles[p].position[1] ÷ PRECISION
            col = particles[p].position[2] ÷ PRECISION
            updateFlow(flow, flow_copy, particle_locs, row, col, cols, 0) # TODO: Implement this.
            particles[p].oldFlow = flow[row, col]
        end
        for p in 1:n_particles
            row = particles[p].position[1] ÷ PRECISION
            col = particles[p].position[2] ÷ PRECISION
            res = particles[p].resistance
            
            back = (particles[p].oldFlow * res ÷ PRECISION) ÷ particle_locs[row, col]
            flow[row, col] -= back
            flow[row - 1, col] += back ÷ 2
            if col > 1
                flow[row - 1, col - 1] += back ÷ 4
            else
                flow[row - 1, col] += back ÷ 4
            end
            if col < cols - 1
                flow[row - 1, col + 1] += back ÷ 4
            else
                flow[row - 1, col] += back ÷ 4
            end
        end
    end
end

function copyToAncillary(rows, cols, flow_copy, flow)
    for r in 1:iter
        if r > rows
            break
        end
        for c in 1:cols
            flow_copy[r, c] = flow[r, c]
        end
    end
end

function propagationStage(iter, threshold, rows, cols, flow, flow_copy, particle_locs)
    # Initialize data to detect max variability.
    if iter % STEPS == 1
        threshold == 0
    end
    # Execute propagation on the wave fronts.
    wave_front = iter % STEPS
    if wave_front == 0
        wave_front = STEPS
    end
    for wave in wave_front:STEPS:rows
        if wave > iter
            continue
        end
        for col in 1:cols
            var = updateFlow(flow, flow_copy, particle_locs, wave, col, cols, 1)
            if var > threshold
                threshold = var
            end
        end
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
    
    # Random sequences initializer.
    random_seq = [parse(Int, ARGS[i]) for i in 13:15]
    seed = rand(random_seq)
    Random.seed!(seed)

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
    # Add number of fixed and moving particles in the bands.
    n_fixed_band = trunc(Int, fixed_band_size * cols * fixed_density)
    n_moving_band = trunc(Int, moving_band_size * cols * moving_density)
    n_particles = n_particles + n_fixed_band + n_moving_band
    particles::Vector{Particle} = []                    # List to store cells information
    
    if rank == 0; printDEBUG(n_fixed_band, n_moving_band, n_particles) end;MPI.Barrier(comm_cart)
    
    # Read fixed particles.
    fixed_particles = (length(ARGS) - 15) ÷ 3
    lowerbound = 1
    for i in lowerbound:fixed_particles
        row = parse(Int, ARGS[13 + i * 3]) * PRECISION
        col = parse(Int, ARGS[14 + i * 3]) * PRECISION
        mass = 0
        resistance = trunc(Int, parse(Float16, ARGS[15 + i * 3]) * PRECISION)
        speed = [0, 0]
        push!(particles, Particle([row, col], mass, resistance, speed, 0))
    end
    
    # Generate fixed particles in band.
    lowerbound = fixed_particles + 1
    for i in lowerbound:n_particles-n_moving_band
        row = trunc(Int, PRECISION * (fixed_band_pos + fixed_band_size * rand(1)[1]))
        col = trunc(Int, PRECISION * cols * rand(1)[1])
        mass = 0
        resistance = trunc(Int, PRECISION * rand(1)[1])
        speed = [0, 0]
        push!(particles, Particle([row, col], mass, resistance, speed, 0))
    end

    # Generate moving particles in band.
    lowerbound = n_particles - n_moving_band + 1
    # if rank == 0; println(lowerbound) end;MPI.Barrier(comm_cart)
    for i in lowerbound:n_particles
        row = trunc(Int, PRECISION * (moving_band_pos + moving_band_size * rand(1)[1]))
        col = trunc(Int, PRECISION * cols * rand(1)[1])
        mass = trunc(Int, PRECISION * (1 + 5 * rand(1)[1]))
        resistance = trunc(Int, PRECISION * rand(1)[1])
        speed = [0, 0]
        push!(particles, Particle([row, col], mass, resistance, speed, 0))
    end

    if rank == 0
        printArguments(rows, cols, max_iter, threshold, fan_pos, fan_size,
            fixed_band_pos, fixed_band_size, fixed_density, 
            moving_band_pos, moving_band_size, moving_density, random_seq, n_particles, particles)
    end
    MPI.Barrier(comm_cart)

    # TODO: Start global timer.

    # Initialization for parallelization.
    flow = zeros(Int, cols, rows)            # Tunnel air flow.
    flow_copy = zeros(Int, cols, rows)       # Tunnel air flow (ancillary, of secondary importance, copy).
    particle_locs = zeros(Int, cols, rows)   # Quickly locate particle positions.

    # Simulation.
    for iter in 1:max_iter
        # Change fan values each STEP iterations.
        updateFan(iter, fan_pos, fan_size, flow)
        
        # Particles movement each STEPS iterations.
        particleMovements(iter, rows, cols, particle_locs, particles, flow)
        
        # Effects due to particles each STEPS iterations.
        particleEffects(iter, cols, n_particles, particles, flow, flow_copy, particle_locs)

        # Copy data in ancillary structure.
        copyToAncillary(rows, cols, flow_copy, flow)

        # Propagation stage.
        propagationStage(iter, threshold, rows, cols, flow, flow_copy, particle_locs)
    end

    # TODO: Stop global timer.

    # TODO: Print results.
end

main()
