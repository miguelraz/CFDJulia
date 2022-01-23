# Text provided under a Creative Commons Attribution license, CC-BY.  All code is made available under the FSF-approved BSD-3 license.  (c) Lorena A. Barba, Gilbert F. Forsyth 2017. Thanks to NSF for support via CAREER award #1149784.
# [@LorenaABarba](https://twitter.com/LorenaABarba)
# Translation to Julia by [@miguelraz_](https://twitter.com/miguelraz_)

# 12 steps to Navier–Stokes
# =====
# ***

# Up to now, all of our work has been in one spatial dimension (Steps [1](./01_Step_1.ipynb) to [4](./05_Step_4.ipynb)). We can learn a lot in just 1D, but let's grow up to flatland: two dimensions.

# In the following exercises, you will extend the first four steps to 2D. To extend the 1D finite-difference formulas to partial derivatives in 2D or 3D, just apply the definition: a partial derivative with respect to $x$ is the variation in the $x$ direction *at constant* $y$.

# In 2D space, a rectangular (uniform) grid is defined by the points with coordinates:

# $$x_i = x_0 +i \Delta x$$

# $$y_i = y_0 +i \Delta y$$

# Now, define $u_{i,j} = u(x_i,y_j)$ and apply the finite-difference formulas on either variable $x,y$ *acting separately* on the $i$ and $j$ indices. All derivatives are based on the 2D Taylor expansion of a mesh point value around $u_{i,j}$.

# Hence, for a first-order partial derivative in the $x$-direction, a finite-difference formula is:

# $$ \frac{\partial u}{\partial x}\biggr\rvert_{i,j} = \frac{u_{i+1,j}-u_{i,j}}{\Delta x}+\mathcal{O}(\Delta x)$$

# and similarly in the $y$ direction. Thus, we can write backward-difference, forward-difference or central-difference formulas for Steps 5 to 12. Let's get started!

# Step 5: 2-D Linear Convection
# ----
# ***

# The PDE governing 2-D Linear Convection is written as

# $$\frac{\partial u}{\partial t}+c\frac{\partial u}{\partial x} + c\frac{\partial u}{\partial y} = 0$$

# This is the exact same form as with 1-D Linear Convection, except that we now have two spatial dimensions to account for as we step forward in time.

# Again, the timestep will be discretized as a forward difference and both spatial steps will be discretized as backward differences.

# With 1-D implementations, we used $i$ subscripts to denote movement in space (e.g. $u_{i}^n-u_{i-1}^n$).  Now that we have two dimensions to account for, we need to add a second subscript, $j$, to account for all the information in the regime.

# Here, we'll again use $i$ as the index for our $x$ values, and we'll add the $j$ subscript to track our $y$ values.

# With that in mind, our discretization of the PDE should be relatively straightforward.

# $$\frac{u_{i,j}^{n+1}-u_{i,j}^n}{\Delta t} + c\frac{u_{i, j}^n-u_{i-1,j}^n}{\Delta x} + c\frac{u_{i,j}^n-u_{i,j-1}^n}{\Delta y}=0$$

# As before, solve for the only unknown:

# $$u_{i,j}^{n+1} = u_{i,j}^n-c \frac{\Delta t}{\Delta x}(u_{i,j}^n-u_{i-1,j}^n)-c \frac{\Delta t}{\Delta y}(u_{i,j}^n-u_{i,j-1}^n)$$

# We will solve this equation with the following initial conditions:

# $$u(x,y) = \begin{cases}
# \begin{matrix}
# 2\ \text{for} & 0.5 \leq x, y \leq 1 \cr
# 1\ \text{for} & \text{everywhere else}\end{matrix}\end{cases}$$

# and boundary conditions:

# $$u = 1\ \text{for } \begin{cases}
# \begin{matrix}
# x =  0,\ 2 \cr
# y =  0,\ 2 \end{matrix}\end{cases}$$


using Plots

### Variable declarations

nₓ = 81
ny = 81 # note we can't call this n\_y because of the Unicode standard, not because of Julia
nₜ = 100
c = 1
Δx = 2 / (nx - 1)
Δy = 2 / (ny - 1)
σ = 0.2
Δt = σ * Δx
rangex = range(start = 0, stop = 2, length = nₓ)
rangey = range(start = 0, stop = 2, length = ny)

### TODO
u = ones(ny, nₓ) ##create a 1xn vector of 1's
uₙ = ones(ny, nₓ)

###Assign initial conditions

##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
# We're going to use an array comprehension for this:
u = [ .5 <= x <= 1 && .5 <= y <= 1 ? 2.0 : 1.0 for x in rangex, y in rangey ]

# Array comprehensions are super useful when initializing data, but might not be as efficient
# They work by putting a `for` expression inside square brackets

firstten = [i for i in 1:10]
squares = [i^2 for i in 1:10]

# You can also filter the elements with a predicate
evensquares = [i^2 for i in 1:10 if iseven(i)]

# and have multiple loops
squaremat = [i^2 for i in 1:4, j in 1:4]

###Plot Initial Condition
### TODO
##the figsize parameter can be used to produce different sized images
fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')                      
X, Y = numpy.meshgrid(x, y)                            
surf = ax.plot_surface(X, Y, u[:], cmap=cm.viridis)

### 3D Plotting Notes

# To plot a projected 3D result, make sure that you have added the Axes3D library.

#     from mpl_toolkits.mplot3d import Axes3D

# The actual plotting commands are a little more involved than with simple 2d plots.

# ```python
# fig = pyplot.figure(figsize=(11, 7), dpi=100)
# ax = fig.gca(projection='3d')
# surf2 = ax.plot_surface(X, Y, u[:])
# ```

# The first line here is initializing a figure window.  The **figsize** and **dpi** commands are optional and simply specify the size and resolution of the figure being produced.  You may omit them, but you will still require the

#     fig = pyplot.figure()

# The next line assigns the plot window the axes label 'ax' and also specifies that it will be a 3d projection plot.  The final line uses the command

#     plot_surface()

# which is equivalent to the regular plot command, but it takes a grid of X and Y values for the data point positions.


# ##### Note


# The `X` and `Y` values that you pass to `plot_surface` are not the 1-D vectors `x` and `y`.  In order to use matplotlibs 3D plotting functions, you need to generate a grid of `x, y` values which correspond to each coordinate in the plotting frame.  This coordinate grid is generated using the numpy function `meshgrid`.

#     X, Y = numpy.meshgrid(x, y)

# ### Iterating in two dimensions

# To evaluate the wave in two dimensions requires the use of several nested for-loops to cover all of the `i`'s and `j`'s.
# ~~~Since Python is not a compiled language there can be noticeable slowdowns in the execution of code with multiple for-loops.~~~
# Since Julia is great at for loops, we'll just write simple for loops
# First try evaluating the 2D convection code and see what results it produces.


# TODO
u = ones((ny, nx))
# u[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2
uₙ = similar(u)

for n in 1:nₜ ##loop across number of time steps
    copyto!(uₙ, u)
    row, col = size(u)
    # TODO C vs FORTRAN order
    for j in 2:row
        for i in 2:col
            u[j, i] = (uₙ[j, i] - (c * Δt / Δx * (uₙ[j, i] - uₙ[j, i - 1])) -
                                  (c * Δt / Δy * (uₙ[j, i] - uₙ[j - 1, i])))
            # Can you see why using begin/end can help to see the rectangular boundary conditions?
            u[begin, :] .= 1.0
            u[end, :]   .= 1.0
            u[:, begin] .= 1.0
            u[:, end]   .= 1.0
        end
    end
end


# TODO
fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
surf2 = ax.plot_surface(X, Y, u[:], cmap=cm.viridis)


# Array Operations
# ----------------

# Here the same 2D convection code is implemented, but instead of using nested for-loops, the same calculations are evaluated using array views.


# TODO timing
u = ones((ny, nx))
u[CartesianIndices((rangex, rangey))] .= 2
uₙ = similar(u)

for n in 1:nₜ ##loop across number of time steps
    copyto!(uₙ, u)
    u[2:end, 2:end] = @views (uₙ[2:end, 2:end] - (c * Δt / Δx * (uₙ[2:end, 2:end] - uₙ[2:end, begin:(end-1)])) -
                              (c * Δt / Δy * (uₙ[2:end, 2:end] - uₙ[begin:(end-1), 2:end])))
    u[begin, :] .= 1.0
    u[end,   :] .= 1.0
    u[:, begin] .= 1.0
    u[:,   end] .= 1.0
end

### Index agnosticism
# The expression `(begin+1):end` represents the same as `2:end`, but only when all the indices are linear. If someone were to shift the indices of `u`
# because it's a special type of matrix (as in OffsetArrays.jl), our code would still work without a hitch or performance hit.
# This is precisely the sort of composability that lets you reuse your code and drag-and-drop algorithms to solve interesting problems
# without having to change your code internals.
# #TODO Unreasonable effectiveness of multiple dispatch.

# TODO
fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
surf2 = ax.plot_surface(X, Y, u[:], cmap=cm.viridis)




## Learn More

The video lesson that walks you through the details for Step 5 (and onwards to Step 8) is **Video Lesson 6** on You Tube:

# TODO
# YouTubeVideo('tUg_dE3NXoY')
# src="https://www.youtube.com/embed/tUg_dE3NXoY"
