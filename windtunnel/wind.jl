using MPI
using Random

PRECISION = 10000
STEPS = 8

# Structure to represent a solid particle in the tunnel surface.
mutable struct Particle
    position::Vector # As [row, column]
    mass::Int
    resistance::Int
    speed::Vector # As [row, column]
    oldFlow::Int
end

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

printDEBUG(n_fixed_band, n_moving_band, n_particles) = println("### printDEBUG ###\n$(length(ARGS)), $(length(ARGS)-15), $((length(ARGS)-15)÷3)\nn_fixed_band=$n_fixed_band\nn_moving_band=$n_moving_band\nn_particles=$n_particles")

globalStartRow(rows, dims, coord) = rows ÷ dims[2] * coord[2]

globalStartColumn(cols, dims, coord) = cols ÷ dims[1] * coord[1]

inCartGrid(gsr, gsc, row, col, own_row, own_col) = (gsr <= row && row < gsr + own_row && gsc <= col && col < gsc + own_col) ? true : false

function readFixedParticles(lb, ub, gsr, gsc, own_row, own_col, particles)
    for i in lb:ub
        row = parse(Int, ARGS[13 + i * 3])
        col = parse(Int, ARGS[14 + i * 3])
        if !inCartGrid(gsr, gsc, row, col, own_row, own_col)
            continue
        end
        mass = 0
        resistance = trunc(Int, parse(Float16, ARGS[15 + i * 3]) * PRECISION)
        speed = [0, 0]
        push!(particles, Particle([row * PRECISION, col * PRECISION], mass, resistance, speed, 0))
    end
end

function generateFixedBand(lb, ub, fixed_band_pos, fixed_band_size, own_col, particles)
    for i in lb:ub
        row = trunc(Int, PRECISION * (fixed_band_pos + fixed_band_size * rand(1)[1]))
        col = trunc(Int, PRECISION * cols * rand(1)[1])
        mass = 0
        resistance = trunc(Int, PRECISION * rand(1)[1])
        speed = [0, 0]
        push!(particles, Particle([row, col], mass, resistance, speed, 0))
    end
end

function generateMovingBand(lb, ub, moving_band_pos, moving_band_size, own_col, particles)
    for i in lb:ub
        row = trunc(Int, PRECISION * (moving_band_pos + moving_band_size * rand(1)[1]))
        col = trunc(Int, PRECISION * cols * rand(1)[1])
        mass = trunc(Int, PRECISION * (1 + 5 * rand(1)[1]))
        resistance = trunc(Int, PRECISION * rand(1)[1])
        speed = [0, 0]
        push!(particles, Particle([row, col], mass, resistance, speed, 0))
    end
end

"""Given a size the function returns the dimensions for 
the cartesian grid as [columns, rows]\n
|0        |1        |
|:--:     |:--:     |
|0 (0, 0) | 1 (1, 0)|

|0        |1        |
|:--:     |:--:     |
|1 (0, 1) | 3 (1, 1)|
|0 (0, 0) | 2 (1, 0)|

|0        |1        |
|:--:     |:--:     |
|3 (0, 3) | 7 (1, 3)|
|2 (0, 2) | 6 (1, 2)|
|1 (0, 1) | 5 (1, 1)|
|0 (0, 0) | 4 (1, 0)|"""
function getDims(size)
    if size == 1
        return [1, 1]
    elseif size == 2
        return [2, 1]
    elseif size == 4
        return [2, 2]
    elseif size == 8
        return [2, 4]
    elseif size == 16
        return [4, 4]
    end
end

"""Update only the first row of the ranks at the top of the topology."""
function updateFan(iter, fan_pos, fan_size, flow, dims, coord, target_row, cols, own_col)
    if !(dims[2] == coord[2] + 1) || !(iter % STEPS == 1); return end

    # block = cols ÷ dims[1] * coord[1]
    gsc = globalStartColumn(cols, dims, coord)
    
    for j in 2:own_col+1
        if (fan_pos <= gsc + j) && (gsc + j <= fan_pos + fan_size)
            phase = iter ÷ STEPS * (π ÷ 4)
            phase_step = π ÷ 2 ÷ fan_size
            pressure_level = 9 + 2 * sin(phase + (j - fan_pos) * phase_step)
            noise = 0.5 - rand(1)[1]
            # Store level in the first row of the ancillary structure.
            flow[target_row, j] = trunc(Int, PRECISION * (pressure_level + noise))
        end
    end
end

# TODO: Check whether the placements of c and r are correct.
function moveParticle(flow, particles, p, rows, cols)
    for i in 1:STEPS
        r = particles[p].position[1] ÷ PRECISION
        c = particles[p].position[2] ÷ PRECISION
        pressure = flow[r - 1, c]
        left = 0; right = 0
        c == 2 ? left = 0 : left = pressure - flow[r - 1, c - 1];
        c == own_col + 1 ? right = 0 : right = pressure - flow[r - 1, c + 1];

        flow_row = trunc(Int, pressure ÷ particles[p].mass * PRECISION)
        flow_col = trunc(Int, (right - left) ÷ particles[p].mass * PRECISION)

        # Speed change.
        particles[p].speed[1] = (particles[p].speed[1] + flow_row) ÷ 2
        particles[p].speed[2] = (particles[p].speed[2] + flow_col) ÷ 2
        
        # Movement.
        particles[p].position[1] = particles[p].position[1] + particles[p].speed[1] ÷ STEPS ÷ 2
        particles[p].position[2] = particles[p].position[2] + particles[p].speed[2] ÷ STEPS ÷ 2
    
        # Control limits.
        if particles[p].position[1] >= PRECISION * own_row # TODO: left here
            particles[p] = PRECISION * rows - 1
        end
        if particles[p].position[1] < 0
            particles[p].position[1] = 0
        end
        if particles[p].position[1] >= PRECISION * cols
            particles[p].position[1] = PRECISION * cols
        end
    end
end

function particleMovements(iter, rows, cols, particle_locs, particles, flow)
    if !(iter % STEPS == 1); return end
    
    # Clean particle positions.
    for r in 2:own_row + 1
        if r > iter
            break
        end
        for c in 2:own_col + 1
            particle_locs[r, c] = 0
        end
    end

    #= 
    TODO: Left here.
    Choose between the following, list of particles to contain 
    entire rows x cols grid or own_row x own_col grid.
    =#
    ub = length(particles)
    for p in 1:ub

        if () && particles[p].mass == 0
            continue
        end
        moveParticle(flow, particles, p, rows, cols)
    end

    # Annotate position.
    for p in 1:n_particles
        row = particles[p].position[1] ÷ PRECISION
        col = particles[p].position[2] ÷ PRECISION + 1
        particle_locs[row, col] += 1
    end
end

function updateFlow(flow, flow_copy, particle_locs, row, col, cols, skip_particles)
    # Skip update in particle positions.
    if skip_particles && (particle_locs[c,r] != 0)
        return 0
    end

    # Update if border left.
    if col == 0
        flow[col, row] = (flow_copy[col, row] + flow_copy[col-1, row] * 2 + flow_copy[col-1, row+1]) ÷ 4
    end
    # Update if border right.
    if col == cols - 1
        flow[col, row] = (flow_copy[col, row] + flow_copy[col-1, row] * 2 + flow_copy[col-1, row-1]) ÷ 4      
    end
    if col > 0 && (col < cols - 1)
        flow[col, row] = (flow_copy[col, row] + flow_copy[col-1, row] * 2 + flow_copy[col-1, row-1] + flow_copy[col-1, row+1]) ÷ 5
    end

    # Return flow variation at this position.
    return abs(flow_copy[col, row] - flow(col, row))
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

Order of arguments: rows, cols, max_iter, threshold, fan_pos, fan_size, fixed_band_pos, fixed_band_size, fixed_density, moving_band_pos, moving_band_size, moving_density, rand_seq_1, rand_seq_2, rand_seq_3, particle_row, particle_col, particle_density, ...
"""
function main()
    # Simulation data.
    rows = parse(Int, ARGS[1])                  # Tunnel rows.
    cols = parse(Int, ARGS[2])                  # Tunnel columns.
    max_iter = parse(Int, ARGS[3])              # Max number of sim (simulation) steps.
    threshold = parse(Float16, ARGS[4])         # Threshold of variability to continue the sim.
    fan_pos = parse(Int, ARGS[5])               # First position of the fan.
    fan_size = parse(Int, ARGS[6])              # Fan size.
    fixed_band_pos = parse(Int, ARGS[7])        # First position of the band where FIXED particles Start.
    fixed_band_size = parse(Int, ARGS[8])       # Size of the band where FIXED particles start.
    fixed_density = parse(Float16, ARGS[9])     # Density of starting FIXED particles.
    moving_band_pos = parse(Int, ARGS[10])      # First position of the band where MOVING particles start.
    moving_band_size = parse(Int, ARGS[11])     # Size of the band where MOVING particles start.
    moving_density = parse(Float16, ARGS[12])   # Density of starting MOVING particles.
    
    # Initialize random sequences.
    random_seq = [parse(Int, ARGS[i]) for i in 13:15]
    seed = rand(random_seq)
    Random.seed!(seed)

    # Initialize MPI.
    MPI.Init()
    comm = MPI.Comm_dup(MPI.COMM_WORLD)
    n_ranks = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)

    # Initialize MPI cartesian grid.
    dims = getDims(n_ranks)
    periodic = map(_ -> false, dims)
    reorder = false
    cart_comm = MPI.Cart_create(comm, dims; periodic, reorder)
    cart_coord = MPI.Cart_coords(cart_comm)

    # Initialize own rows and cols. dims = [c, r].
    own_col = cols ÷ dims[1]
    own_row = rows ÷ dims[2]

    # Read particles from the arguments.
    n_particles = (length(ARGS) - 15) ÷ 3
    # Add number of fixed and moving particles in the bands.
    n_fixed_band = trunc(Int, fixed_band_size * cols * fixed_density)
    n_moving_band = trunc(Int, moving_band_size * cols * moving_density)
    n_particles = n_particles + n_fixed_band + n_moving_band
    # List to store cells information.
    particles::Vector{Particle} = []
    
    # if rank == 0; printDEBUG(n_fixed_band, n_moving_band, n_particles) end;MPI.Barrier(cart_comm)
    
    # Read fixed particles only if it is in corresponding cartesian grid.
    fixed_particles = (length(ARGS) - 15) ÷ 3
    gsr = globalStartRow(rows, dims, cart_coord)
    gsc = globalStartColumn(cols, dims, cart_coord)
    readFixedParticles(1, fixed_particles, gsr, gsc, own_row, own_col, particles)
    
    # Generate fixed particles in band.
    generateFixedBand(fixed_particles + 1, n_particles-n_moving_band, fixed_band_pos, fixed_band_size, own_col, particles)

    # Generate moving particles in band.
    generateMovingBand(n_particles - n_moving_band + 1, n_particles, moving_band_pos, moving_band_size, own_col, particles)

    if rank == 0; printArguments(rows, cols, max_iter, threshold, fan_pos, fan_size, fixed_band_pos, fixed_band_size, fixed_density, moving_band_pos, moving_band_size, moving_density, random_seq, n_particles, particles) end
    
    # TODO: Start global timer.
    MPI.Barrier(cart_comm)

    # Initialization for parallelization.
    flow = zeros(Int, own_row + 2, own_col + 2)            # With ghostcell. Tunnel air flow.
    flow_copy = zeros(Int, own_row + 2, own_col + 2)       # With ghostcell. Tunnel air flow (ancillary, of secondary importance, copy).
    particle_locs = zeros(Int, own_row + 2, own_col + 2)   # With ghostcell. Quickly locate particle positions.

    # Initialize ghostcells with value -1.
    flow[1, :] = flow[end, :] .= -1; flow[:, 1] = flow[:, end] .= -1
    flow_copy[1, :] = flow_copy[end, :] .= -1; flow_copy[:, 1] = flow_copy[:, end] .= -1
    particle_locs[1, :] = particle_locs[end, :] .= -1; particle_locs[:, 1] = particle_locs[:, end] .= -1

    # Simulation.
    for iter in 1:max_iter
        # Change fan values each STEP iterations.
        updateFan(iter, fan_pos + 1, fan_size, flow, dims, cart_coord, 1, cols, own_col)

        # # Particles movement each STEPS iterations.
        # particleMovements(iter, rows, cols, particle_locs, particles, flow)
        
        # # Effects due to particles each STEPS iterations.
        # particleEffects(iter, cols, n_particles, particles, flow, flow_copy, particle_locs)

        # # Copy data in ancillary structure.
        # copyToAncillary(rows, cols, flow_copy, flow)

        # # Propagation stage.
        # propagationStage(iter, threshold, rows, cols, flow, flow_copy, particle_locs)
    end

    # TODO: Stop global timer.
    MPI.Barrier(cart_comm)

    # TODO: Print results.
end

# Comment out call to main() whenever running tests.jl.
main()
