using Random
using BenchmarkTools

PRECISION = 10000
STEPS = 8


mutable struct Particle
    position::Vector
    mass::Int
    resistance::Int
    speed::Vector
    old_flow::Int
end


function read_fixed_particles(lb, ub, particles)
    for i in lb:ub
        r = parse(Int, ARGS[13 + i * 3]) * PRECISION
        c = parse(Int, ARGS[14 + i * 3]) * PRECISION
        mass = 0
        resistance = trunc(Int, parse(Float16, ARGS[15 + i * 3]) * PRECISION)
        speed = [0, 0]
        push!(particles, Particle([r, c], mass, resistance, speed, 0))
    end
end


function generate_fixed_particles(lb, ub, particles, fbp, fbs, tc)
    for i in lb:ub
       r = trunc(Int, PRECISION * (fbp + fbs * rand(1)[1]))
       c = trunc(Int, PRECISION * tc * rand(1)[1])
       mass = 0
       resistance = trunc(Int, PRECISION * rand(1)[1])
       speed = [0, 0]
       push!(particles, Particle([r, c], mass, resistance, speed, 0))
    end
end


function generate_moving_particles(lb, ub, particles, mbp, mbs, tc)
    for i in lb:ub
        r = trunc(Int, PRECISION * (mbp + mbs * rand(1)[1]))
        c = trunc(Int, PRECISION * tc * rand(1)[1])
        mass = trunc(Int, PRECISION * (1 + 5 * rand(1)[1]))
        resistance = trunc(Int, PRECISION * rand(1)[1])
        speed = [0, 0]
        push!(particles, Particle([r, c], mass, resistance, speed, 0)) 
    end
end


function update_fan(iter, fp , fs, flow)
    if iter % STEPS != 1
        return
    end

    # println(" update_fan start.")
    for c in fp:fp + fs - 1
        phase = iter ÷ STEPS * (π ÷ 4)
        phase_step = π ÷ 2 ÷ fs
        pressure = 9 + 2 * sin(phase + (c - fp) * phase_step)
        noise = 0.5 - rand(1)[1]
        flow[1, c] = trunc(Int, PRECISION * (pressure + noise))
    end

    # println(" update_fan exit.")
end


function p_move(flow, p, tr, tc)
    # println("   p_move start.")
    for step in 0:STEPS
        # println("     step = $step ")
        r = p.position[1] ÷ PRECISION
        c = p.position[2] ÷ PRECISION == 0 ? 1 : p.position[2] ÷ PRECISION
        pressure = flow[r - 1, c]

        # c == 1 ? left = 0 : left = pressure - flow[r - 1, c - 1];
        # c == tc ? right = 0 : right = pressure - flow[r - 1, c + 1];

        if c == 1
            left = 0
        else
            left = pressure - flow[r - 1, c - 1]
        end

        if c == tc
            right = 0
        else
            right = pressure - flow[r - 1, c + 1]
        end

        flow_r = trunc(Int, pressure ÷ p.mass * PRECISION)
        flow_c = trunc(Int, (right - left) ÷ p.mass * PRECISION)

        # Speed change.
        p.speed[1] = (p.speed[1] + flow_r) ÷ 2
        p.speed[2] = (p.speed[2] + flow_c) ÷ 2

        # Movement.
        p.position[1] = p.position[1] + p.speed[1] ÷ STEPS ÷ 2
        p.position[2] = p.position[2] + p.speed[2] ÷ STEPS ÷ 2

        # Boundaries.
        if p.position[1] >= PRECISION * tr
            p.position[1] = PRECISION * tr
        end
        if p.position[2] <= 1
            p.position[2] = 1
        end
        if p.position[2] >= PRECISION * tc
            p.position[2] = PRECISION * tc
        end
    end
    # println("   p_move exit.")
end


function p_movements(iter, tr, tc, pl, particles, flow)
    if iter % STEPS != 1
        return
    end

    # println(" p_movements start.")

    # Clean particle locs.
    for r in 1:tr
        if r > iter
            break
        end
        for c in 1:tc
            pl[r, c] = 0
        end
    end

    # Move particles.
    # println(" particles = $particles ")
    for p in particles
        if p.mass == 0
            continue
        end
        p_move(flow, p, tr, tc)
    end

    # Annotate particle locations.
    for p in particles
        r = p.position[1] ÷ PRECISION
        c = p.position[2] ÷ PRECISION == 0 ? 1 : p.position[2] ÷ PRECISION
        pl[r, c] += 1
    end
    # println(" p_movements exit.")
end


function update_flow(flow, flow_copy, pl, r, c, tc, skip)
    # println(" update_flow start.")
    if skip && pl[r, c] != 0 || r == 1
        return 0
    end
    # Update border left.
    if c == 1
        flow[r, c] = (flow_copy[r, c] + flow_copy[r - 1, c] * 2 + flow_copy[r - 1, c + 1]) ÷ 4
    end
    # Update border right.
    if c == tc
        flow[r, c] = (flow_copy[r, c] + flow_copy[r - 1, c] * 2 + flow_copy[r - 1, c - 1]) ÷ 4
    end
    # Update central part.
    if 1 < c && c < tc
        flow[r, c] = (flow_copy[r, c] + flow_copy[r - 1, c] * 2 + flow_copy[r - 1, c - 1] + flow_copy[r - 1, c + 1]) ÷ 5
    end
    # println(" update_flow exit.")
    return abs(flow_copy[r, c] - flow[r, c])
end


function p_effects(iter, particles, flow, flow_copy, p_locs, tc)
    if iter % STEPS != 1
        return
    end

    # println("  p_effects start.")

    for p in particles
        r = p.position[1] ÷ PRECISION
        c = p.position[2] ÷ PRECISION == 0 ? 1 : p.position[2] ÷ PRECISION
        update_flow(flow, flow_copy, p_locs, r, c, tc, false)
        p.old_flow = flow[r, c]
    end

    for p in particles
        r = p.position[1] ÷ PRECISION
        c = p.position[2] ÷ PRECISION == 0 ? 1 : p.position[2] ÷ PRECISION
        back = trunc(Int, (p.old_flow * p.resistance ÷ PRECISION) ÷ p_locs[r, c])
        flow[r, c] -= back
        flow[r - 1, c] += back ÷ 2

        c > 1 ? flow[r - 1, c - 1] += back ÷ 4 : flow[r - 1, c] += back ÷ 4;
        c < tc ? flow[r - 1, c + 1] += back ÷ 4 : flow[r - 1, c] += back ÷ 4;
    end
    # println("  p_effects exit.")
end


function p_wavefront(iter, max_var, tr, tc, flow, flow_copy, p_locs)
    # println(" p_wavefront start.")
    wavefront = iter % STEPS
    if wavefront == 1
        max_var = 0
    elseif wavefront == 0
        wavefront = STEPS
    end

    wave = wavefront
    while wave < tr
        # println("   wave = $wave ")
        if wave > iter
            break
        end

        for c in 1:tc
            var = update_flow(flow, flow_copy, p_locs, wave, c, tc, true)
            if var > max_var
                max_var = var
            end
        end

        wave += STEPS
    end
    # println(" p_wavefront exit.")
end


function simulate_windtunnel(start_iter, max_iter, max_var, threshold, fan_pos, fan_size, tr, tc, flow, flow_copy, p_locs, particles)
    # println("START simulation.")
    for iter in start_iter:max_iter
        # println(" iter = $iter ")
        if max_var <= threshold
            break
        end

        update_fan(iter, fan_pos, fan_size, flow)

        p_movements(iter, tr, tc, p_locs, particles, flow)

        p_effects(iter, particles, flow, flow_copy, p_locs, tc)

        copy!(flow_copy, flow)

        p_wavefront(iter, max_var, tr, tc, flow, flow_copy, p_locs)
    end
    # println("EXIT simulation.")
end


function print_status(tr, tc, flow, pl)
    # println(" print_status start. ")
    result = Matrix{String}(undef, tr, tc)

    for r in 1:tr, c in 1:tc
        if flow[r, c] >= 10 * PRECISION
            symbol = '*'
        elseif flow[r, c] >= 1 * PRECISION
            symbol = '0' + flow[r, c] ÷ PRECISION;
        elseif flow[r, c] >= 0.5 * PRECISION
            symbol = '+'
        elseif flow[r, c] > 0
            symbol = '.'
        else
            symbol = ' '
        end

        pl[r, c] > 0 ? result[r, c] = "[$symbol]" : result[r, c] = string(symbol)
    end

    display(result)

    res = ""
    for r in 1:tr
        for c in 1:tc
            length(result[r, c]) > 1 ? res *= " $(result[r,c]) " : res *= "  $(result[r,c])  "
        end
        res *= "\n"
    end

    write("restult.txt", res)

    # println(" print_status exit. ")
end


function main()
    # println("start main")

    # Simulation data.
    tunnel_rows = parse(Int, ARGS[1])
    tunnel_cols = parse(Int, ARGS[2])
    max_iter = parse(Int, ARGS[3])              # Max number of simulation steps.
    threshold = parse(Float16, ARGS[4])         # Threshold of variability to continue the simulation.
    fan_pos = parse(Int, ARGS[5])               # Fan position in the tunnel.
    fan_size = parse(Int, ARGS[6])              # Fan size in the tunnel.

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

    remaining_args = (length(ARGS) - 15) ÷ 3
    n_p_f_b = trunc(Int, fixed_band_size * tunnel_cols * fixed_density)
    n_p_m_b = trunc(Int, moving_band_size * tunnel_cols * moving_density)
    n_particles = remaining_args + n_p_f_b + n_p_m_b
    particles::Vector{Particle} = []

    read_fixed_particles(1, remaining_args, particles)
    ub = n_particles - n_p_m_b - remaining_args
    generate_fixed_particles(1, ub, particles, fixed_band_pos, fixed_band_size, tunnel_cols)
    ub = n_particles - n_p_f_b - remaining_args
    generate_moving_particles(1, ub, particles, moving_band_pos, moving_band_size, tunnel_cols)

    max_var = typemax(Int64)

    flow = zeros(Int, tunnel_rows, tunnel_cols)
    flow_copy = copy(flow)
    p_locs = copy(flow)

    start_iter = 0
    @time simulate_windtunnel(start_iter, max_iter, max_var, threshold, fan_pos, fan_size, tunnel_rows, tunnel_cols, flow, flow_copy, p_locs, particles)

    print_status(tunnel_rows, tunnel_cols, flow, p_locs)
end

main()