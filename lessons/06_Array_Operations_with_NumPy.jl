# Text provided under a Creative Commons Attribution license, CC-BY.  All code is made available under the FSF-approved BSD-3 license.  (c) Lorena A. Barba, Gilbert F. Forsyth 2017. Thanks to NSF for support via CAREER award #1149784.
# [@LorenaABarba](https://twitter.com/LorenaABarba)
# Translation to Julia by [@miguelraz_](https://twitter.com/miguelraz_)

## 12 steps to Navier–Stokes
# =====
# ***

# This lesson complements the first interactive module of the online [CFD Python](https://github.com/barbagroup/CFDPython) class, by Prof. Lorena A. Barba, called **12 Steps to Navier–Stokes.** It was written with BU graduate student Gilbert Forsyth.

## Array Operations with NumPy
# ----------------

# For more computationally intensive programs, the use of built-in Numpy functions can provide an  increase in execution speed many-times over.  As a simple example, consider the following equation:

# $$u^{n+1}_i = u^n_i-u^n_{i-1}$$

# Now, given a vector $u^n = [0, 1, 2, 3, 4, 5]\ \ $   we can calculate the values of $u^{n+1}$ by iterating over the values of $u^n$ with a for loop.


u = [0, 1, 2, 3, 4, 5]

for i in 1:length(u)
    println(u[i] - u[i-1])
end


# This is the expected result and the execution time was nearly instantaneous.
# We can calculate each operation with one function:

diff(u)

# What this command says is subtract the 1st, 2nd, 3rd, 4th and 5th elements of $u$ from the 2nd, 3rd, 4th, 5th and 6th elements of $u$.

### Speed Increases
# In NumPy, we would have to learn a bucketful about how to restructure our program to exploit "vectorization".
# This unfortunately makes code very brittle and diminishes readability, but it's require by NumPy for speed.
# We have no such obstacles in Julia - if you write simple for loops, they will start out being "fast enough",
# and then we can move on to talk about worthy optimizations later.
# Note: There


nₓ = 81
ny = 81 # note we can't call this n\_y because of the Unicode standard, not because of Julia
nₜ = 100
c = 1
Δx = 2 / (nx - 1)
Δy = 2 / (ny - 1)
σ = 0.2
Δt = σ * Δx

x = range(start = 0, stop = 2, length = nₓ)
y = range(start = 0, stop = 2, length = ny)

### TODO
u = ones(ny, nₓ) ##create a 1xn vector of 1's
uₙ = ones(ny, nₓ)

###Assign initial conditions
#We can try the following to start the initial conditions:
# TODO - cartesian indexes
rangey = (.5/Δy):(1/Δy)
rangex = (.5/Δx):(1/Δx)
for (_, j) in rangey
    for (_, i) in rangex
        u[i, j] = 2.0
    end
end
u

# But we could also try doing this with `CartesianIndexes` - they will subset a rectangular region of our array (for any dimension!)
# and loop efficiently through it. You can find more about them if you press `?CartesianIndex` in the REPL
#u[int(.5 / dy): int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2

u[CartesianIndices((rangex,rangey))] .= 2

# With our initial conditions all set up, let's first try running our original nested loop code, making use of the
# BenchmarkTools.jl `@btime` macro, which will help us evaluate the performance of our code.

# **Note**: The `@btime` macro will run the code several times and then give an average execution time as a result.
# If you have any figures being plotted within a code cell where you run `@btime`,
# it will plot those figures repeatedly which can be a bit messy and slow.

# The execution times below will vary from machine to machine.  Don't expect your times to match these times,
# but you _should_ expect to see the same general trend in decreasing execution time as we switch to array operations.

using BenchmarkTools

# TODO
@btime u[int(.5 / dy): int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2 :setup=(u = ones(ny, nx))


@btime begin
    for n in 1:nₜ ##loop across number of time steps
    uₙ = copy(u)
    row, col = size(u)
    for j in 1:row
        for i in 1:col
            u[j, i] = (uₙ[j, i] - (c * Δt / Δx *
                                  (uₙ[j, i] - uₙ[j, i - 1])) -
                                  (c * Δt / Δy *
                                   (uₙ[j, i] - uₙ[j - 1, i])))
            u[1, :] = 1.0
            u[end, :] = 1.0
            u[:, 1] = 1.0
            u[:, end] = 1.0
        end
    end
    end
end :setup=(u = ones(ny, nx))

## Python comparison
# With the "raw" Python code above, the mean execution time achieved was 3.07 seconds (on a MacBook Pro Mid 2012).
# Keep in mind that with these three nested loops, that the statements inside the **j** loop are being evaluated more than 650,000 times.
# Let's compare that with the performance of the same code implemented with array operations:
# I've left the vectorized version so that people can compare the code and **appreciate** how substantial the difference from Julia is.


# ```python
# %%timeit
# u = numpy.ones((ny, nx))
# u[int(.5 / dy): int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2

# for n in range(nt + 1): ##loop across number of time steps
#     un = u.copy()
#     u[1:, 1:] = (un[1:, 1:] - (c * dt / dx * (un[1:, 1:] - un[1:, 0:-1])) -
#                               (c * dt / dy * (un[1:, 1:] - un[0:-1, 1:])))
#     u[0, :] = 1
#     u[-1, :] = 1
#     u[:, 0] = 1
#     u[:, -1] = 1
# ```

#     7.38 ms ± 105 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)


# As you can see, the speed increase is substantial in Python.  The same calculation goes from 3.07 seconds to 7.38 milliseconds.
# 3 seconds isn't a huge amount of time to wait, but these speed gains will increase exponentially
# with the size and complexity of the problem being evaluated.

## Julia timings:
# We're going to be doing some mild modifications to the code, which will involve
# - avoid global variables (even though Julia 1.8 will not penalize this), as it's good practice
# - diminishing allocations (because [views in Julia copy by default](TODO perf manual and views))
# - putting code inside of a function
# There's very little point to write efficient array code if it's not within a function in Julia.
# This is because Julia needs to JIT compile your code to get all the speedups.
# The JIT compiler does not need you to use type hints or annotations to work properly, it's only necessary that
# [TODO type stability]() function output types be predictable on input types only.
#
function time1(u, params)
    u = ones(ny, nₓ)
    (;Δt, Δx, c, nₜ) = params

    for n in 1:nₜ ##loop across number of time steps
    uₙ = copy(u)
            @views u[:,:] .= (uₙ[2:end, 2:end] - (c * Δt / Δx *
                                  (uₙ[2:end, 2:end] - uₙ[2:end, 1:end])) -
                                  (c * Δt / Δy *
                                   (uₙ[2:end, 2:end] - uₙ[1:end, 2:end])))
            @views u[1, :] .= 1.0
            @views u[end, :] .= 1.0
            @views u[:, 1] .= 1.0
            @views u[:, end] .= 1.0
        end
    end
    end
end
params = (;Δt, Δx, c, nₜ)
@btime time1(u, params)
