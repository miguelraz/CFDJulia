# Text provided under a Creative Commons Attribution license, CC-BY.  All code is made available under the FSF-approved BSD-3 license.  (c) Lorena A. Barba, Gilbert F. Forsyth 2017. Thanks to NSF for support via CAREER award #1149784.
# [@LorenaABarba](https://twitter.com/LorenaABarba)
# Translation to Julia by [@miguelraz_](https://twitter.com/miguelraz_)

# 12 steps to Navier–Stokes
# =====
# ***

# This lesson complements the first interactive module of the online [CFD Python](https://github.com/barbagroup/CFDPython) class,
# by Prof. Lorena A. Barba, called **12 Steps to Navier–Stokes.** The interactive module starts with
# simple exercises in 1D that at first use little of the power of Python.
# We now present some new ways of doing the same things that are more efficient and produce prettier code, but in Julia 🚀 ,
# which is hosted at the [CFD Julia](TODO LINK) repo

# This lesson was written with BU graduate student Gilbert Forsyth.


# Defining Functions in Julia
# ----

# In steps 1 through 8, we wrote Python code that is meant to run from top to bottom.
# We were able to reuse code (to great effect!) by copying and pasting, to incrementally build a solver for the Burgers' equation.
# But moving forward there are more efficient ways to write our Python codes.
# In this lesson, we are going to introduce *function definitions*,
# which will allow us more flexibility in reusing and also in organizing our code.

# We'll begin with a trivial example: a function which adds two numbers.

# To create a function in Julia, we start with the following function block:

function simpleadd(a,b)
    return a + b
end


# The `return` statement tells Julia what data to return in response to being called.
# In Julia, if you omit the `return` inside the body of a funciton, the last evaluated expression is implicitly returned.
# Now we can try calling our `simpleadd` function:

simpleadd(3, 4)

# For simple functions, we can use the inline function syntax:
simpleadd(a,b) = a + b

# Where both are equivalent.

# Of course, there can be much more happening between the `function` line and the `return` line.
# In this way, one can build code in a *modular* way.
# Let's try a function which returns the `n`-th number in the Fibonacci sequence.

function fibonacci(n)
    a, b = 0, 1
    for i in 1:n
        a, b = b, a + b
    end
    return a
end

fibonacci(7)

# Once defined, the function `fibonacci` can be called like any of the built-in Python functions that we've already used.
# For exmaple, we might want to print out the Fibonacci sequence up through the `n`-th value:


for n in 1:10
    println(fibonacci(n))
end


# We will use the capacity of defining our own functions in Julia to help us build code that is easier to reuse, easier to maintain, easier to share!

##### Exercise

# (Pending.)

# Learn more
# -----
# ***

# Remember our short detour on using [array operations with Julia](./07_Step_5.jl)?

# Well, there are a few more ways to make your scientific codes in Julia run faster.
# The original tutorial said...
# > We recommend the article on the Technical Discovery blog about [Speeding Up Python](http://technicaldiscovery.blogspot.com/2011/06/speeding-up-python-numpy-cython-and.html) (June 20, 2011), which talks about NumPy, Cython and Weave. It uses as example the Laplace equation (which we will solve in [Step 9](./12_Step_9.ipynb)) and makes neat use of defined functions.
#  We *encourage* instead that you read [the performance tips section of the Julia manual](https://docs.julialang.org/en/v1/manual/performance-tips/).

# From the original tutorial, it says...
# > But a recent new way to get fast Python codes is [Numba](http://numba.pydata.org).
# Which is actually a good comparison to how Julia's JIT operates internally, but unfortunately breaks down in many cases
# [that we are interested in](http://www.stochasticlifestyle.com/why-numba-and-cython-are-not-substitutes-for-julia/)
