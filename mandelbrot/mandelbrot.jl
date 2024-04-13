#=
Exercise 3 from: https://www.francescverdugo.com/XM_40017/dev/julia_basics/

Function mandel estimates if a given point (x,y) in the complex plane belongs to the Mandelbrot set.

If the value of mandel is less than max_iters, the point is provably outside the Mandelbrot set.
If mandel is equal to max_iters, then the point is provably inside the set.
The larger max_iters, the better the quality of the estimate (the nicer will be your plot).

Plot the value of function mandel for each pixel in a 2D grid of the box.

(−1.7, 0.7) × (−1.2, 1.2)

Use a grid resolution of at least 1000 points in each direction and max_iters at least 10.
You can increase these values to get nicer plots.
To plot the values use function heatmap from the Julia package GLMakie.
Use LinRange to divide the horizontal and vertical axes into pixels.
See the documentation of these functions for help.
GLMakie is a GPU-accelerated plotting back-end for Julia.
It is a large package and it can take some time to install and to generate the first plot.
Be patient.
=#

using GLMakie

MAX_ITERS = 1000
GRID_RESOLUTION = 10000
xmin = -1.7
xmax = 0.7
ymin = -1.2
ymax = 1.2

function mandel(x, y, max_iters)
    z = Complex(x, y)
    c = z
    threshold = 2
    for n in 1:max_iters
        if abs(z) > threshold
            return n - 1
        end
        z = z^2 + c
    end
    max_iters
end

function plotMandel()
    xAxis = LinRange(xmin, xmax, GRID_RESOLUTION)
    yAxis = LinRange(ymin, ymax, GRID_RESOLUTION)
    
    data = [mandel(x, y, MAX_ITERS) for x in xAxis, y in yAxis]
    p = heatmap(data)
    save("mandelbrot.png", p)
end

plotMandel()