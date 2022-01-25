# WIP - This repo is incomplete and needs a brave contributor like you!

> Please see the issues to see how you can contribute.

# CFD Julia

> Please cite as: Barba, Lorena A., and Forsyth, Gilbert F. (2018). CFD Python: the 12 steps to Navier-Stokes equations. _Journal of Open Source Education_, **1**(9), 21, https://doi.org/10.21105/jose.00021

# [![DOI](https://jose.theoj.org/papers/10.21105/jose.00021/status.svg)](https://doi.org/10.21105/jose.00021)

**CFD Python**, a.k.a. the **12 steps to Navier-Stokes**, is a practical module for learning the foundations of Computational Fluid Dynamics (CFD) by coding solutions to the basic partial differential equations that describe the physics of fluid flow.

(**CFD Julia**, however, is a Julia port of the same tutorial so as to not run into the numerous headachaes that Python brings to high performance computing. It is a port done by Miguel Raz Guzmán, and the comments are surround by `()` or wstart with `NB`. )

The module was part of a course taught by [Prof. Lorena Barba](http://lorenabarba.com) between 2009 and 2013 in the Mechanical Engineering department at Boston University (Prof. Barba since moved to the George Washington University).

The module assumes only basic programming knowledge (in any language) and some background in partial differential equations and fluid mechanics. The "steps" were inspired by ideas of Dr. Rio Yokota, who was a post-doc in Prof. Barba's lab until 2011, and the lessons were refined by Prof. Barba and her students over several semesters teaching the CFD course. 
We wrote this set of Jupyter notebooks in 2013 to teach an intensive two-day course in Mendoza, Argentina.
(Miguel rewrote this tutorial to have a side by side comparison between Python and Julia for STEM beginners.)

Guiding students through these steps (without skipping any!), they learn many valuable lessons. The incremental nature of the exercises means they get a sense of achievement at the end of each assignment, and they feel they are learning with low effort. As they progress, they naturally practice code re-use and they incrementally learn programming and plotting techniques. As they analyze their results, they learn about numerical diffusion, accuracy and convergence. 
In about four weeks of a regularly scheduled course, they become moderately proficient programmers and are motivated to start discussing more theoretical matters.

## How to use this module

In a regular-session university course, students can complete the **CFD Python** lessons in 4 to 5 weeks. 
As an intensive tutorial, the module can be completed in two or three full days, depending on the learner's prior experience. 
The lessons can also be used for self study. 
In all cases, learners should follow along the worked examples in each lesson by re-typing the code in a fresh Jupyter notebook, maybe taking original notes as they try things out. 
(You may run the vanilla Julia files or the interactive Pluto.jl notebooks, take your pic.)

Lessons
-------
Steps 1–4 are in one spatial dimension. Steps 5–10 are in two dimensions (2D). Steps 11–12 solve the Navier-Stokes equation in 2D. Three "bonus" notebooks cover the CFL condition for numerical stability, array operations with NumPy, and defining functions in Python.

> FIX numerical libraries
> FIX links
* [Quick Julia Intro]() — For Python novices, this lesson introduces the numerical libraries (NumPy and Matplotlib), Python variables, use of whitespace, and slicing arrays.
* [Step 1]() — Linear convection with a step-function initial condition (IC) and appropriate boundary conditions (BCs).
* [Step 2]() — With the same IC/BCs, _nonlinear_ convection.
* [CFL Condition]() — Exploring numerical stability and the Courant-Friedrichs-Lewy (CFL) condition.
* [Step 3]() — With the same IC/BCs, _diffusion_ only.
* [Step 4]() — Burgers’ equation, with a saw-tooth IC and periodic BCs (with an introduction to Sympy).
* [Array Operations with Julia]()
* [Step 5]() — Linear convection in 2D with a square-function IC and appropriate BCs.
* [Step 6]() — With the same IC/BCs, _nonlinear_ convection in 2D.
* [Step 7]() — With the same IC/BCs, _diffusion_ in 2D.
* [Step 8]() — Burgers’ equation in 2D
* [Defining Functions in Julia]() - How to define functions in Julia.
* [Step 9]() — Laplace equation with zero IC and both Neumann and Dirichlet BCs.
* [Step 10]() — Poisson equation in 2D.
* [Step 11]() — Solves the Navier-Stokes equation for 2D cavity flow.
* [Step 12]() — Solves the Navier-Stokes equation for 2D channel flow.



## Dependencies

#### TODO Fix dependencies
To use these lessons, you need Python 3, and the standard stack of scientific Python: NumPy, Matplotlib, SciPy, Sympy. And of course, you need [Jupyter](http://jupyter.org)—an interactive computational environment that runs on a web browser.
To use these lessons, you need Julia 1.6 or above and all the packages you need can be installed by typing
> ] instantiate .

This mini-course is built as a set of [Pluto notebooks]() containing the written materials and worked-out solutions on Julia code. To work with the material, we recommend that you start each lesson with a fresh new notebook, and follow along, typing each line of code (don't copy-and-paste!), and exploring by changing parameters and seeing what happens. 


#### Installing via Anaconda

NOPE! Don't even worry about it. Go take a nap or pet a dog and think about all the problems you don't have.

### Running the notebooks

They are normal Julia files! Or you can also do

```julia
using Pluto
Pluto.run_notebook()
```
And run them from there.

This will start up a Pluto session in your browser!

## How to contribute to CFD Julia

We accept contributions via pull request—in fact, several users have already submitted pull requests making corrections or small improvements. You can also open an issue if you find a bug, or have a suggestion. 

## Copyright and License

(c) 2017 Lorena A. Barba, Gilbert F. Forsyth. All content is under Creative Commons Attribution [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/legalcode.txt), and all [code is under BSD-3 clause](https://github.com/engineersCode/EngComp/blob/master/LICENSE) (previously under MIT, and changed on March 8, 2018). 

We are happy if you re-use the content in any way!

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) [![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

