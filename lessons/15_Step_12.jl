# Text provided under a Creative Commons Attribution license, CC-BY.  All code is made available under the FSF-approved BSD-3 license.  (c) Lorena A. Barba, Gilbert F. Forsyth 2017. Thanks to NSF for support via CAREER award #1149784.
# [@LorenaABarba](https://twitter.com/LorenaABarba)
# Translation to Julia by [@miguelraz_](https://twitter.com/miguelraz_)

# 12 steps to Navier–Stokes
# =====
# ***

# Did you make it this far? This is the last step! How long did it take you to write your own Navier–Stokes solver in Python following this interactive module? Let us know!

# Step 12: Channel Flow with Navier–Stokes
# ----
# ***

# The only difference between this final step and Step 11 is that we are going to add a source term to the $u$-momentum equation, to mimic the effect of a pressure-driven channel flow. Here are our modified Navier–Stokes equations:

# $$\frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=-\frac{1}{\rho}\frac{\partial p}{\partial x}+\nu\left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 u}{\partial y^2}\right)+F$$

# $$\frac{\partial v}{\partial t}+u\frac{\partial v}{\partial x}+v\frac{\partial v}{\partial y}=-\frac{1}{\rho}\frac{\partial p}{\partial y}+\nu\left(\frac{\partial^2 v}{\partial x^2}+\frac{\partial^2 v}{\partial y^2}\right)$$

# $$\frac{\partial^2 p}{\partial x^2}+\frac{\partial^2 p}{\partial y^2}=-\rho\left(\frac{\partial u}{\partial x}\frac{\partial u}{\partial x}+2\frac{\partial u}{\partial y}\frac{\partial v}{\partial x}+\frac{\partial v}{\partial y}\frac{\partial v}{\partial y}\right)
# $$

# ### Discretized equations

# With patience and care, we write the discretized form of the equations. It is highly recommended that you write these in your own hand, mentally following each term as you write it.

# The $u$-momentum equation:

# $$
# \begin{split}
# & \frac{u_{i,j}^{n+1}-u_{i,j}^{n}}{\Delta t}+u_{i,j}^{n}\frac{u_{i,j}^{n}-u_{i-1,j}^{n}}{\Delta x}+v_{i,j}^{n}\frac{u_{i,j}^{n}-u_{i,j-1}^{n}}{\Delta y} = \\
# & \qquad -\frac{1}{\rho}\frac{p_{i+1,j}^{n}-p_{i-1,j}^{n}}{2\Delta x} \\
# & \qquad +\nu\left(\frac{u_{i+1,j}^{n}-2u_{i,j}^{n}+u_{i-1,j}^{n}}{\Delta x^2}+\frac{u_{i,j+1}^{n}-2u_{i,j}^{n}+u_{i,j-1}^{n}}{\Delta y^2}\right)+F_{i,j}
# \end{split}
# $$

# The $v$-momentum equation:

# $$
# \begin{split}
# & \frac{v_{i,j}^{n+1}-v_{i,j}^{n}}{\Delta t}+u_{i,j}^{n}\frac{v_{i,j}^{n}-v_{i-1,j}^{n}}{\Delta x}+v_{i,j}^{n}\frac{v_{i,j}^{n}-v_{i,j-1}^{n}}{\Delta y} = \\
# & \qquad -\frac{1}{\rho}\frac{p_{i,j+1}^{n}-p_{i,j-1}^{n}}{2\Delta y} \\
# & \qquad +\nu\left(\frac{v_{i+1,j}^{n}-2v_{i,j}^{n}+v_{i-1,j}^{n}}{\Delta x^2}+\frac{v_{i,j+1}^{n}-2v_{i,j}^{n}+v_{i,j-1}^{n}}{\Delta y^2}\right)
# \end{split}
# $$

# And the pressure equation:

# $$
# \begin{split}
# & \frac{p_{i+1,j}^{n}-2p_{i,j}^{n}+p_{i-1,j}^{n}}{\Delta x^2} + \frac{p_{i,j+1}^{n}-2p_{i,j}^{n}+p_{i,j-1}^{n}}{\Delta y^2} = \\
# & \qquad \rho\left[\frac{1}{\Delta t}\left(\frac{u_{i+1,j}-u_{i-1,j}}{2\Delta x}+\frac{v_{i,j+1}-v_{i,j-1}}{2\Delta y}\right) - \frac{u_{i+1,j}-u_{i-1,j}}{2\Delta x}\frac{u_{i+1,j}-u_{i-1,j}}{2\Delta x} - 2\frac{u_{i,j+1}-u_{i,j-1}}{2\Delta y}\frac{v_{i+1,j}-v_{i-1,j}}{2\Delta x} - \frac{v_{i,j+1}-v_{i,j-1}}{2\Delta y}\frac{v_{i,j+1}-v_{i,j-1}}{2\Delta y}\right]
# \end{split}
# $$

# As always, we need to re-arrange these equations to the form we need in the code to make the iterations proceed.

# For the $u$- and $v$ momentum equations, we isolate the velocity at time step `n+1`:

# $$
# \begin{split}
# u_{i,j}^{n+1} = u_{i,j}^{n} & - u_{i,j}^{n} \frac{\Delta t}{\Delta x} \left(u_{i,j}^{n}-u_{i-1,j}^{n}\right) - v_{i,j}^{n} \frac{\Delta t}{\Delta y} \left(u_{i,j}^{n}-u_{i,j-1}^{n}\right) \\
# & - \frac{\Delta t}{\rho 2\Delta x} \left(p_{i+1,j}^{n}-p_{i-1,j}^{n}\right) \\
# & + \nu\left[\frac{\Delta t}{\Delta x^2} \left(u_{i+1,j}^{n}-2u_{i,j}^{n}+u_{i-1,j}^{n}\right) + \frac{\Delta t}{\Delta y^2} \left(u_{i,j+1}^{n}-2u_{i,j}^{n}+u_{i,j-1}^{n}\right)\right] \\
# & + \Delta t F
# \end{split}
# $$

# $$
# \begin{split}
# v_{i,j}^{n+1} = v_{i,j}^{n} & - u_{i,j}^{n} \frac{\Delta t}{\Delta x} \left(v_{i,j}^{n}-v_{i-1,j}^{n}\right) - v_{i,j}^{n} \frac{\Delta t}{\Delta y} \left(v_{i,j}^{n}-v_{i,j-1}^{n}\right) \\
# & - \frac{\Delta t}{\rho 2\Delta y} \left(p_{i,j+1}^{n}-p_{i,j-1}^{n}\right) \\
# & + \nu\left[\frac{\Delta t}{\Delta x^2} \left(v_{i+1,j}^{n}-2v_{i,j}^{n}+v_{i-1,j}^{n}\right) + \frac{\Delta t}{\Delta y^2} \left(v_{i,j+1}^{n}-2v_{i,j}^{n}+v_{i,j-1}^{n}\right)\right]
# \end{split}
# $$

# And for the pressure equation, we isolate the term $p_{i,j}^n$ to iterate in pseudo-time:

# $$
# \begin{split}
# p_{i,j}^{n} = & \frac{\left(p_{i+1,j}^{n}+p_{i-1,j}^{n}\right) \Delta y^2 + \left(p_{i,j+1}^{n}+p_{i,j-1}^{n}\right) \Delta x^2}{2(\Delta x^2+\Delta y^2)} \\
# & -\frac{\rho\Delta x^2\Delta y^2}{2\left(\Delta x^2+\Delta y^2\right)} \\
# & \times \left[\frac{1}{\Delta t} \left(\frac{u_{i+1,j}-u_{i-1,j}}{2\Delta x} + \frac{v_{i,j+1}-v_{i,j-1}}{2\Delta y}\right) - \frac{u_{i+1,j}-u_{i-1,j}}{2\Delta x}\frac{u_{i+1,j}-u_{i-1,j}}{2\Delta x} - 2\frac{u_{i,j+1}-u_{i,j-1}}{2\Delta y}\frac{v_{i+1,j}-v_{i-1,j}}{2\Delta x} - \frac{v_{i,j+1}-v_{i,j-1}}{2\Delta y}\frac{v_{i,j+1}-v_{i,j-1}}{2\Delta y}\right]
# \end{split}
# $$

# The initial condition is $u, v, p=0$ everywhere, and at the boundary conditions are:

# $u, v, p$ are periodic on $x=0,2$

# $u, v =0$ at $y =0,2$

# $\frac{\partial p}{\partial y}=0$ at $y =0,2$

# $F=1$ everywhere.

# Let's begin by importing our usual run of libraries:

using Plots

# In step 11, we isolated a portion of our transposed equation to make it easier to parse and we're going to do the same thing here.
# One thing to note is that we have periodic boundary conditions throughout this grid,
# so we need to explicitly calculate the values at the leading and trailing edge of our `u` vector.


function build_up_b(rho, dt, dx, dy, u, v)
    b = zeros(u)
    b[2:end, 1:-1] = (rho * (1 / dt * ((u[2:end, 2:] - u[2:end, 0:-2]) / (2 * dx) +
                                      (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -
                            ((u[2:end, 2:] - u[2:end, 0:-2]) / (2 * dx))^2 -
                            2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) *
                                 (v[2:end, 2:] - v[2:end, 0:-2]) / (2 * dx))-
                            ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))^2))
    
    # Periodic BC Pressure @ x = 2
    b[2:end, -1] = (rho * (1 / dt * ((u[2:end, 0] - u[2:end,-2]) / (2 * dx) +
                                    (v[2:, -1] - v[0:-2, -1]) / (2 * dy)) -
                          ((u[2:end, 0] - u[2:end, -2]) / (2 * dx))^2 -
                          2 * ((u[2:, -1] - u[0:-2, -1]) / (2 * dy) *
                               (v[2:end, 0] - v[2:end, -2]) / (2 * dx)) -
                          ((v[2:, -1] - v[0:-2, -1]) / (2 * dy))^2))

    # Periodic BC Pressure @ x = 0
    b[2:end, 0] = (rho * (1 / dt * ((u[2:end, 1] - u[2:end, -1]) / (2 * dx) +
                                   (v[2:, 0] - v[0:-2, 0]) / (2 * dy)) -
                         ((u[2:end, 1] - u[2:end, -1]) / (2 * dx))^2 -
                         2 * ((u[2:, 0] - u[0:-2, 0]) / (2 * dy) *
                              (v[2:end, 1] - v[2:end, -1]) / (2 * dx))-
                         ((v[2:, 0] - v[0:-2, 0]) / (2 * dy))^2))
    
    return b
end

# We'll also define a Pressure Poisson iterative function, again like we did in Step 11.
# Once more, note that we have to include the periodic boundary conditions at the leading and trailing edge.
# We also have to specify the boundary conditions at the top and bottom of our grid.


function pressure_poisson_periodic(p, dx, dy)
    pn = empty(p)
    
    for q in range(nit):
        pn = p.copy()
        p[2:end, 1:-1] = (((pn[2:end, 2:] + pn[2:end, 0:-2]) * dy^2 +
                          (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx^2) /
                         (2 * (dx^2 + dy^2)) -
                         dx^2 * dy^2 / (2 * (dx^2 + dy^2)) * b[2:end, 1:-1])

        # Periodic BC Pressure @ x = 2
        p[2:end, -1] = (((pn[2:end, 0] + pn[2:end, -2])* dy^2 +
                        (pn[2:, -1] + pn[0:-2, -1]) * dx^2) /
                       (2 * (dx^2 + dy^2)) -
                       dx^2 * dy^2 / (2 * (dx^2 + dy^2)) * b[2:end, -1])

        # Periodic BC Pressure @ x = 0
        p[2:end, 0] = (((pn[2:end, 1] + pn[2:end, -1])* dy^2 +
                       (pn[2:, 0] + pn[0:-2, 0]) * dx^2) /
                      (2 * (dx^2 + dy^2)) -
                      dx^2 * dy^2 / (2 * (dx^2 + dy^2)) * b[2:end, 0])
        
        # Wall boundary conditions, pressure
        p[-1, :] =p[-2, :]  # dp/dy = 0 at y = 2
        p[0, :] = p[1, :]  # dp/dy = 0 at y = 0
        end
    
    return p
end

Now we have our familiar list of variables and initial conditions to declare before we start.


# Variable declarations
nx = 41
ny = 41
nt = 10
nit = 50 
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
x = LinRange(0, 2, nx)
y = LinRange(0, 2, ny)
X, Y = numpy.meshgrid(x, y)


##physical variables
rho = 1
nu = .1
F = 1
dt = .01

#initial conditions
u = zeros((ny, nx))
un = zeros((ny, nx))

v = zeros((ny, nx))
vn = zeros((ny, nx))

p = ones((ny, nx))
pn = ones((ny, nx))

b = zeros((ny, nx))

# For the meat of our computation, we're going to reach back to a trick we used in Step 9 for Laplace's Equation.
# We're interested in what our grid will look like once we've reached a near-steady state.
# We can either specify a number of timesteps `nt` and increment it until we're satisfied with the results,
# or we can tell our code to run until the difference between two consecutive iterations is very small.

# We also have to manage **8** separate boundary conditions for each iteration.
# The code below writes each of them out explicitly.  If you're interested in a challenge,
# you can try to write a function which can handle some or all of these boundary conditions.
# TODO Fix data structure recommendation
# If you're interested in tackling that, you should probably read up on Julia [dictionaries](http://docs.python.org/2/tutorial/datastructures.html#dictionaries).


udiff = 1
stepcount = 0

while udiff > .001
    un = copy(u)
    vn = copy(v)

    b = build_up_b(rho, dt, dx, dy, u, v)
    p = pressure_poisson_periodic(p, dx, dy)

    u[2:end, 2:end] = (un[2:end, 2:end] -
                     un[2:end, 2:end] * dt / dx *
                    (un[2:end, 2:end] - un[2:end, begin:end-1]) -
                     vn[2:end, 2:end] * dt / dy *
                    (un[2:end, 2:end] - un[begin:end-1, 2:end]) -
                     dt / (2 * rho * dx) * 
                    (p[2:end, 3:end] - p[2:end, begin:end-1]) +
                     nu * (dt / dx^2 *
                    (un[2:end, 3:end] - 2 * un[2:end, 2:end] + un[2:end, begin:end-1]) +
                     dt / dy^2 *
                    (un[3:end, 2:end] - 2 * un[2:end, 2:end] + un[begin:end-1, 2:end])) +
                     F * dt)

    v[2:end, 2:end] = (vn[2:end, 2:end] -
                     un[2:end, 2:end] * dt / dx *
                    (vn[2:end, 2:end] - vn[2:end, begin:end-1]) -
                     vn[2:end, 2:end] * dt / dy *
                    (vn[2:end, 2:end] - vn[begin:end-1, 2:end]) -
                     dt / (2 * rho * dy) * 
                    (p[3:end, 2:end] - p[begin:end-1, 2:end]) +
                     nu * (dt / dx^2 *
                    (vn[2:end, 3:end] - 2 * vn[2:end, 2:end] + vn[2:end, begin:end-1]) +
                     dt / dy^2 *
                    (vn[3:end, 2:end] - 2 * vn[2:end, 2:end] + vn[begin:end-1, 2:end])))

    # Periodic BC u @ x = 2     
    u[2:end, end] = (un[2:end, end] - un[2:end, end] * dt / dx *
                  (un[2:end, end] - un[2:end, end-1]) -
                   vn[2:end, end] * dt / dy *
                  (un[2:end, end] - un[begin:end-1, end]) -
                   dt / (2 * rho * dx) *
                  (p[2:end, begin] - p[2:end, end-1]) +
                   nu * (dt / dx^2 *
                  (un[2:end, begin] - 2 * un[2:end,-1] + un[2:end, end-1]) +
                   dt / dy^2 *
                  (un[3:end, end] - 2 * un[2:end, end] + un[begin:end-1, end])) + F * dt)

    # Periodic BC u @ x = 0
    u[2:end, begin] = (un[2:end, begin] - un[2:end, begin] * dt / dx *
                 (un[2:end, begin] - un[2:end, end]) -
                  vn[2:end, begin] * dt / dy *
                 (un[2:end, begin] - un[begin:end-1, begin]) -
                  dt / (2 * rho * dx) * 
                 (p[2:end, 2] - p[2:end, end]) +
                  nu * (dt / dx^2 *
                 (un[2:end, 2] - 2 * un[2:end, begin] + un[2:end, end]) +
                  dt / dy^2 *
                 (un[3:end, begin] - 2 * un[2:end, begin] + un[begin:end-1, begin])) + F * dt)

    # Periodic BC v @ x = 2
    v[2:end, end] = (vn[2:end, end] - un[2:end, end] * dt / dx *
                  (vn[2:end, end] - vn[2:end, end-1]) -
                   vn[2:end, end] * dt / dy *
                  (vn[2:end, end] - vn[begin:end-1, end]) -
                   dt / (2 * rho * dy) * 
                  (p[3:end, end] - p[begin:end-1, end]) +
                   nu * (dt / dx^2 *
                  (vn[2:end, begin] - 2 * vn[2:end, end] + vn[2:end, end-1]) +
                   dt / dy^2 *
                  (vn[3:end, end] - 2 * vn[2:end, end] + vn[begin:end-1, end])))

    # Periodic BC v @ x = 0
    v[2:end, begin] = (vn[2:end, begin] - un[2:end, begin] * dt / dx *
                 (vn[2:end, begin] - vn[2:end, end]) -
                  vn[2:end, begin] * dt / dy *
                 (vn[2:end, begin] - vn[begin:end-1, begin]) -
                  dt / (2 * rho * dy) * 
                 (p[3:end, begin] - p[begin:end-1, begin]) +
                  nu * (dt / dx^2 *
                 (vn[2:end, 2] - 2 * vn[2:end, begin] + vn[2:end, end]) +
                  dt / dy^2 *
                 (vn[3:end, begin] - 2 * vn[2:end, begin] + vn[begin:end-1, begin])))


    # Wall BC: u,v = 0 @ y = 0,2
    u[begin, :] .= 0
    u[end,   :] .= 0
    v[begin, :] .= 0
    v[end,   :] .= 0

    udiff = (sum(u) - sum(un)) / sum(u)
    stepcount += 1
end

# You can see that we've also included a variable `stepcount` to see how many iterations our loop went through before our stop condition was met.


print(stepcount)

# If you want to see how the number of iterations increases as our `udiff` condition gets smaller and smaller,
# try defining a function to perform the `while` loop written above that takes an input
# `udiff` and outputs the number of iterations that the function runs.

# For now, let's look at our results.  We've used the quiver function
# # to look at the cavity flow results and it works well for channel flow, too.


#TODO PLOT
fig = pyplot.figure(figsize = (11,7), dpi=100)
pyplot.quiver(X[::3, ::3], Y[::3, ::3], u[::3, ::3], v[::3, ::3]);




# The structures in the `quiver` command that look like `[::3, ::3]` are useful when dealing with large amounts of data that you want to visualize.  The one used above tells `matplotlib` to only plot every 3rd data point.  If we leave it out, you can see that the results can appear a little crowded.


# TODO PLOT
fig = pyplot.figure(figsize = (11,7), dpi=100)
pyplot.quiver(X, Y, u, v);


## Learn more
# ***

##### What is the meaning of the $F$ term?

# Step 12 is an exercise demonstrating the problem of flow in a channel or pipe.
# If you recall from your fluid mechanics class, a specified
# pressure gradient is what drives Poisseulle flow.

# Recall the $x$-momentum equation:

# $$\frac{\partial u}{\partial t}+u \cdot \nabla u = -\frac{\partial p}{\partial x}+\nu \nabla^2 u$$

# What we actually do in Step 12 is split the pressure into steady and unsteady components $p=P+p'$. The applied steady pressure gradient is the constant $-\frac{\partial P}{\partial x}=F$ (interpreted as a source term), and the unsteady component is $\frac{\partial p'}{\partial x}$. So the pressure that we solve for in Step 12 is actually $p'$, which for a steady flow is in fact equal to zero everywhere.

# <b>Why did we do this?</b>

# Note that we use periodic boundary conditions for this flow. For a flow with a constant pressure gradient, the value of pressure on the left edge of the domain must be different from the pressure at the right edge. So we cannot apply periodic boundary conditions on the pressure directly. It is easier to fix the gradient and then solve for the perturbations in pressure.

# <b>Shouldn't we always expect a uniform/constant $p'$ then?</b>

# That's true only in the case of steady laminar flows. At high Reynolds numbers, flows in channels can become turbulent, and we will see unsteady fluctuations in the pressure, which will result in non-zero values for $p'$.

# In step 12, note that the pressure field itself is not constant, but it's the pressure perturbation field that is. The pressure field varies linearly along the channel with slope equal to the pressure gradient. Also, for incompressible flows, the absolute value of the pressure is inconsequential.


# ##### And explore more CFD materials online

# The interactive module **12 steps to Navier–Stokes** is one of several components
# of the Computational Fluid Dynamics class taught by Prof. Lorena A. Barba in Boston University between 2009 and 2013.

# For a sample of what the othe components of this class are, you can explore
# the **Resources** section of the Spring 2013 version of [the course's Piazza site](https://piazza.com/bu/spring2013/me702/resources).

# ***


