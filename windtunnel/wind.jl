using MPI

PRECISION = 10000
STEPS = 8

# Structure to represent a solid particle in the tunnel surface.
struct Particle
    position::Vector
    mass::Int
    resistance::Int
    speed::Vector
    oldFlow::Int
end

function main()
    max_iter = 1
end

# TODO: Function to get wall time.

# TODO: Function to simplify matrix/2D array access.

# TODO: Update flow in a matrix position.

# TODO: Move particle function.

# TODO: Print the current state of the simulation.

# TODO: Print usage line in stderr.

# TODO: Main program.

# TODO: Start global timer

# TODO: Parallel initialization.

# TODO: Simulation.

# TODO: Particle movement each STEPS iterations.

# TODO: Effects due to particles each STEPS iterations.

# TODO: Propagation stage.

# TODO: MPI: Fill result arrays used for later output.

# TODO: Stop global timer.

# TODO: Output for the leaderboard.
