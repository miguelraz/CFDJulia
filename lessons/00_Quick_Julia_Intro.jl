# # Julia Crash Course


# Hello! This is a quick intro to programming in Julia to help you hit the ground running with the _12 Steps to Navier–Stokes_.  

# There are two ways to enjoy these lessons with Julia:

# 1. You can download and install a Julia distribution on your computer. Our recommendation is via [julialang.org](julialang.org). Open them up in your favorite IDE and run them one at a time.

# 2. You can run  Julia in a Pluto notebook with `using Pluto; Pluto.run()`.

# In either case, you will probably want to download a copy of this notebook, or the whole AeroPython (?) collection. We recommend that you then follow along each lesson, experimenting with the code in the notebooks, or typing the code into a separate Julia interactive session.

#If you decided to work on your local Julia installation, you will have to navigate in the terminal to the folder that contains the `.jl` files. Then, to launch the notebook server, just type:
#`julia> using Pluto; Pluto.run()`

#You will get a new browser window or tab with a list of the notebooks available in that folder. Click on one and start working!

## Libraries

#Julia is a high-level open-source language.  But the _Julia world_ is inhabited by many packages or libraries that provide useful things like gpu kernels, plotting functions, and much more. We can import libraries of functions to expand the capabilities of Python in our programs.  

#OK! We'll start by importing a few libraries to help us out. First: our favorite library is **NumPy**, ... Just kidding! We don't need NumPy. We will need a 2D plotting library which we will to plot our results.
#The following code will be at the top of most of your programs, so execute this cell first:


#```julia
## <-- comments in Julia are denoted by the pound sign, like this one

#using LinearAlgebra                # we import the array library
#using Plots    # import plotting library
#```

#We are importing one library named `LinearAlgebra` and we are importing a package called `Plots`.
#To use a function belonging to one of these libraries, we have to tell Julia where to look for it. For that, each function name is written following the library name, with a dot in between.
#So if we want to use the unexported functions in `LinearAlgebra`, we have to put a name in front of it:
#
#

myarray =  [ 1 2; 3 4 ]
myarray
LinearAlgebra.norm(myarray)
5.47722

# In this tutorial, it's not too important to know where functions are defined since `using Pkg` calls them all into the global namespace to be used, but if
# you were developing your own code to share with others some additional details about namespacing could save you some headaches.

## Variables

# Julia doesn't require explicitly declared variable types like C and other languages.


a = 5        #a is an integer 5
b = "five"   #b is a string of the word "five'
c = 5.0      #c is a floating point 5  


type(a)
type(b)
type(c)

# Note that if you divide an integer by an integer that yields a float. If you want to use integer division use the `div` operator.

## Whitespace in Julia
#
# Doesn't exist!
#
## Slicing Arrays

#In Julia, you can look at portions of arrays in the same way as in `Matlab`, with a few extra tricks thrown in.  Let's take an array of values from 1 to 5.


myvals = [1, 2, 3, 4, 5]

#Julia uses a **one-based index**, so let's look at the first and last element in the array `myvals`


myvals[1], myvals[5]


#Note here, the slice is inclusive on the front end the back end.

## Assigning Array Variables

# One of the strange little quirks/features in Python that often confuses people comes up when assigning and comparing arrays of values.
# Here is a quick example.  Let's start by defining a 1-D array called $a$:

a = collect(range(start = 1.0, stop = 5.0, length = 5))

#OK, so we have an array $a$, with the values 1 through 5.  I want to make a copy of that array, called $b$, so I'll try the following:

b = a

# Great.  So $a$ has the values 1 through 5 and now so does $b$.
# Now that I have a backup of $a$, I can change its values without worrying about losing data (or so I may think!).


a[2] = 17
a

#Here, the 2nd element of $a$ has been changed to 17.  Now let's check on $b$.
b

# And that's how things go wrong!  When you use a statement like $a = b$, rather than copying all the values of $a$ into a new array called $b$,
# Julia just creates an alias (or a pointer) called $b$ and tells it to route us to $a$.
# So if we change a value in $a$ then $b$ will reflect that change (technically, this is called *assignment by reference*).
# If you want to make a true copy of the array, you have to tell Julia to copy every element of $a$ into a new array.  Let's call it $c$.


c = copy(a)

# Now, we can try again to change a value in $a$ and see if the changes are also seen in $c$.


a[2] = 3
a
c

# OK, it worked!  If the difference between `a = b` and `a = b.copy()` is unclear, you should read through this again.  This issue will come back to haunt you otherwise.

## Learn More

# There are a lot of resources online to learn more about using Julia and other libraries.
# Just for kicks, here we use Pluto's feature for embedding videos to point you to a short video on YouTube on using NumPy arrays.


## TODO
#from IPython.display import YouTubeVideo
# a short video about using NumPy arrays, from Enthought
#YouTubeVideo('vWkb7VahaXQ')
