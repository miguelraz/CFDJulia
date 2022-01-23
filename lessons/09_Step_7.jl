# Text provided under a Creative Commons Attribution license, CC-BY.  All code is made available under the FSF-approved BSD-3 license.  (c) Lorena A. Barba, Gilbert F. Forsyth 2017. Thanks to NSF for support via CAREER award #1149784.
# [@LorenaABarba](https://twitter.com/LorenaABarba)
# Translation to Julia by [@miguelraz_](https://twitter.com/miguelraz_)

# 12 steps to Navier–Stokes
# =====
# ***

# You see where this is going ... we'll do 2D diffusion now and next we will combine steps 6 and 7 to solve Burgers' equation. So make sure your previous steps work well before continuing.

# Step 7: 2D Diffusion
# ----
# ***

# And here is the 2D-diffusion equation:

# $$\frac{\partial u}{\partial t} = \nu \frac{\partial ^2 u}{\partial x^2} + \nu \frac{\partial ^2 u}{\partial y^2}$$

# You will recall that we came up with a method for discretizing second order derivatives in Step 3, when investigating 1-D diffusion.  We are going to use the same scheme here, with our forward difference in time and two second-order derivatives.

# $$\frac{u_{i,j}^{n+1} - u_{i,j}^n}{\Delta t} = \nu \frac{u_{i+1,j}^n - 2 u_{i,j}^n + u_{i-1,j}^n}{\Delta x^2} + \nu \frac{u_{i,j+1}^n-2 u_{i,j}^n + u_{i,j-1}^n}{\Delta y^2}$$

# Once again, we reorganize the discretized equation and solve for $u_{i,j}^{n+1}$

# $$
# \begin{split}
# u_{i,j}^{n+1} = u_{i,j}^n &+ \frac{\nu \Delta t}{\Delta x^2}(u_{i+1,j}^n - 2 u_{i,j}^n + u_{i-1,j}^n) \\
# &+ \frac{\nu \Delta t}{\Delta y^2}(u_{i,j+1}^n-2 u_{i,j}^n + u_{i,j-1}^n)
# \end{split}
# $$

using Plots

# ## Variable declarations
nₓ = 31
ny = 31
nₜ = 17
nᵤ = .05
c = 1
Δx = 2 / (nₓ - 1)
Δy = 2 / (ny - 1)
σ = 0.25
Δt = σ * Δx

rangex = range(start = 0, stop = 2, length = nₓ)
rangey = range(start = 0, stop = 2, length = ny)
indices = CartesianIndices((rangex, rangey))

u = ones(ny, nₓ) ##create a 1xn vector of 1's
uₙ = similar(u)

# ## Assign initial conditions
##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
u[indices] .= 2.0
##set hat function I.C. : v(.5<=x<=1 && .5<=y<=1 ) is 2
v[indices] .= 2.0

# TODO PLOT
fig = pyplot.figure()
ax = fig.gca(projection='3d')
X, Y = numpy.meshgrid(x, y)
surf = ax.plot_surface(X, Y, u, rstride=1, cstride=1, cmap=cm.viridis,
        linewidth=0, antialiased=False)

ax.set_xlim(0, 2)
ax.set_ylim(0, 2)
ax.set_zlim(1, 2.5)

ax.set_xlabel('$x$')
ax.set_ylabel('$y$');



# $$
# \begin{split}
# u_{i,j}^{n+1} = u_{i,j}^n &+ \frac{\nu \Delta t}{\Delta x^2}(u_{i+1,j}^n - 2 u_{i,j}^n + u_{i-1,j}^n) \\
# &+ \frac{\nu \Delta t}{\Delta y^2}(u_{i,j+1}^n-2 u_{i,j}^n + u_{i,j-1}^n)
# \end{split}
# $$


#Run through nt timesteps
# TODO wrap in a let block maybe?
function diffuse(nt)
    rangex = range(start = 0, stop = 2, length = nₓ)
    rangey = range(start = 0, stop = 2, length = ny)
    indices = CartesianIndices((rangex, rangey))
    u[indices] = 2.0
    un = similar(u)

    for n in 1:nt
        copyto!(un, u)
        # TODO check indexing
        u[2:end, 2:end] = (un[2:end,2:end] +
                        nu * Δt / Δx^2 *
                        (un[2:end, 3:end] - 2 * un[2:end, 2:end] + un[2:end, 1:end-1]) +
                        nu * Δt / Δy^2 *
                        (un[3:end, 2:end-1] - 2 * un[2:end, 2:end] + un[begin:end-1, 2:end]))
        u[begin, :] = 1.0
        u[end, :] = 1.0
        u[:, begin] = 1.0
        u[:, end] = 1.0
    end
    u
end



# TODO PLOT
    fig = pyplot.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, u[:], rstride=1, cstride=1, cmap=cm.viridis,
        linewidth=0, antialiased=True)
    ax.set_zlim(1, 2.5)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$');
    



## TODO Plot these 3
diffuse(10)

diffuse(14)

diffuse(50)

## Learn More

The video lesson that walks you through the details for Steps 5 to 8 is **Video Lesson 6** on You Tube:

# TODO YOUTUBE
# src="https://www.youtube.com/embed/tUg_dE3NXoY"
