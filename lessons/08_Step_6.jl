# Text provided under a Creative Commons Attribution license, CC-BY.  All code is made available under the FSF-approved BSD-3 license.  (c) Lorena A. Barba, Gilbert F. Forsyth 2017. Thanks to NSF for support via CAREER award #1149784.
# [@LorenaABarba](https://twitter.com/LorenaABarba)
# Translation to Julia by [@miguelraz_](https://twitter.com/miguelraz_)

# 12 steps to Navier–Stokes
# =====
# ***

# You should have completed your own code for [Step 5](./07_Step_5.ipynb) before continuing to this lesson. As with Steps 1 to 4, we will build incrementally, so it's important to complete the previous step!

# We continue ...

# Step 6: 2-D Convection
# ----
# ***

# Now we solve 2D Convection, represented by the pair of coupled partial differential equations below:

# $$\frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} = 0$$

# $$\frac{\partial v}{\partial t} + u \frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} = 0$$

# Discretizing these equations using the methods we've applied previously yields:

# $$\frac{u_{i,j}^{n+1}-u_{i,j}^n}{\Delta t} + u_{i,j}^n \frac{u_{i,j}^n-u_{i-1,j}^n}{\Delta x} + v_{i,j}^n \frac{u_{i,j}^n-u_{i,j-1}^n}{\Delta y} = 0$$

# $$\frac{v_{i,j}^{n+1}-v_{i,j}^n}{\Delta t} + u_{i,j}^n \frac{v_{i,j}^n-v_{i-1,j}^n}{\Delta x} + v_{i,j}^n \frac{v_{i,j}^n-v_{i,j-1}^n}{\Delta y} = 0$$

# Rearranging both equations, we solve for $u_{i,j}^{n+1}$ and $v_{i,j}^{n+1}$, respectively.  Note that these equations are also coupled.

# $$u_{i,j}^{n+1} = u_{i,j}^n - u_{i,j} \frac{\Delta t}{\Delta x} (u_{i,j}^n-u_{i-1,j}^n) - v_{i,j}^n \frac{\Delta t}{\Delta y} (u_{i,j}^n-u_{i,j-1}^n)$$

# $$v_{i,j}^{n+1} = v_{i,j}^n - u_{i,j} \frac{\Delta t}{\Delta x} (v_{i,j}^n-v_{i-1,j}^n) - v_{i,j}^n \frac{\Delta t}{\Delta y} (v_{i,j}^n-v_{i,j-1}^n)$$

# ### Initial Conditions

# The initial conditions are the same that we used for 1D convection, applied in both the x and y directions.

# $$u,\ v\ = \begin{cases}\begin{matrix}
# 2 & \text{for } x,y \in (0.5, 1)\times(0.5,1) \cr
# 1 & \text{everywhere else}
# \end{matrix}\end{cases}$$

# ### Boundary Conditions

# The boundary conditions hold u and v equal to 1 along the boundaries of the grid
# .

# $$u = 1,\ v = 1 \text{ for } \begin{cases} \begin{matrix}x=0,2\cr y=0,2 \end{matrix}\end{cases}$$

using Plots


### Variable declarations

nₓ = 101
ny = 101
nₜ = 80
c = 1
Δx = 2 / (nₓ - 1)
Δy = 2 / (ny - 1)
σ = 0.2
Δt = σ * Δx

rangex = range(start = 0, stop = 2, length = nₓ)
rangey = range(start = 0, stop = 2, length = ny)
indices = CartesianIndices((rangex, rangey))

uₙ = ones(ny, nₓ) ##create a 1xn vector of 1's
vₙ = ones(ny, nₓ)
u  = ones(ny, nₓ)
v  = ones(ny, nₓ)

###Assign initial conditions
##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
u[indices] .= 2.0
##set hat function I.C. : v(.5<=x<=1 && .5<=y<=1 ) is 2
v[indices] .= 2.0


# TODO PLOT
# Why is [meshgrid](https://groups.google.com/g/julia-users/c/83Pfg9HGhGQ/m/9G_0wi-GBQAJ?pli=1) not an efficient idea?
# Summary: Just use a 2D array comprehension
fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = numpy.meshgrid(x, y)

ax.plot_surface(X, Y, u, cmap=cm.viridis, rstride=2, cstride=2)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$');



# TODO BIKESHED Syntax
# Possible OffsetArrays.jl example?
for n in  1:nₜ
    un = copy(u)
    vn = copy(v)
    @views u[begin+1:end, begin+1:end] = (un[begin+1:end, begin+1:end] -
                 (un[begin+1:end, begin+1:end] * c * Δt / Δx * (un[begin+1:end, begin+1:end] - un[begin+1:end, begin:end-1])) -
                  vn[begin+1:end, begin+1:end] * c * Δt / Δy * (un[begin+1:end, begin+1:end] - un[:-1, begin+1:end]))
    @views v[begin+1:end, begin+1:end] = (vn[begin+1):end, begin+1:end] -
                 (un[begin+1:end, begin+1:end] * c * Δt / Δx * (vn[begin+1:end, begin+1:end] - vn[begin+1:end, begin:end-1])) -
                  vn[begin+1:end, begin+1:end] * c * Δt / Δy * (vn[begin+1:end, begin+1:end] - vn[begin:end-1, begin+1:end]))
    
    u[begin, :] = 1
    u[end, :] = 1
    u[:, begin] = 1
    u[:, end] = 1
    
    v[begin, :] = 1
    v[end, :] = 1
    v[:, begin] = 1
    v[:, end] = 1
end


# TODO PLOT
fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = numpy.meshgrid(x, y)

ax.plot_surface(X, Y, u, cmap=cm.viridis, rstride=2, cstride=2)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$');


# TODO PLOT
fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = numpy.meshgrid(x, y)
ax.plot_surface(X, Y, v, cmap=cm.viridis, rstride=2, cstride=2)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$');


## Learn More

# The video lesson that walks you through the details for Steps 5 to 8 is **Video Lesson 6** on You Tube:


    # TODO YOUTUBE
    # src="https://www.youtube.com/embed/tUg_dE3NXoY"
