# Python Crash Course

Hello! This is a quick intro to programming in Python to help you hit the ground running with the _12 Steps to Navier–Stokes_.  

There are two ways to enjoy these lessons with Python:

1. You can download and install a Python distribution on your computer. One option is the free [Anaconda Scientific Python](https://store.continuum.io/cshop/anaconda/) distribution. Another is [Canopy](https://www.enthought.com/products/canopy/academic/), which is free for academic use.  Our recommendation is Anaconda.

2. You can run Python in the cloud using [Wakari](https://wakari.io/) web-based data analysis, for which you need to create a free account. (No software installation required!)

In either case, you will probably want to download a copy of this notebook, or the whole AeroPython collection. We recommend that you then follow along each lesson, experimenting with the code in the notebooks, or typing the code into a separate Python interactive session.

If you decided to work on your local Python installation, you will have to navigate in the terminal to the folder that contains the .ipynb files. Then, to launch the notebook server, just type:
ipython notebook

You will get a new browser window or tab with a list of the notebooks available in that folder. Click on one and start working!

## Libraries

Python is a high-level open-source language.  But the _Python world_ is inhabited by many packages or libraries that provide useful things like array operations, plotting functions, and much more. We can import libraries of functions to expand the capabilities of Python in our programs.  

OK! We'll start by importing a few libraries to help us out. First: our favorite library is **NumPy**, providing a bunch of useful array operations (similar to MATLAB). We will use it a lot! The second library we need is **Matplotlib**, a 2D plotting library which we will use to plot our results.
The following code will be at the top of most of your programs, so execute this cell first:


```python
# <-- comments in python are denoted by the pound sign, like this one

import numpy                 # we import the array library
from matplotlib import pyplot    # import plotting library
```

We are importing one library named `numpy` and we are importing a module called `pyplot` of a big library called `matplotlib`.
To use a function belonging to one of these libraries, we have to tell Python where to look for it. For that, each function name is written following the library name, with a dot in between.
So if we want to use the NumPy function [linspace()](http://docs.scipy.org/doc/numpy/reference/generated/numpy.linspace.html), which creates an array with equally spaced numbers between a start and end, we call it by writing:


```python
myarray = numpy.linspace(0, 5, 10)
myarray
```




    array([ 0.        ,  0.55555556,  1.11111111,  1.66666667,  2.22222222,
            2.77777778,  3.33333333,  3.88888889,  4.44444444,  5.        ])



If we don't preface the `linspace()` function with `numpy`, Python will throw an error.


```python
myarray = linspace(0, 5, 10)
```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-4-ed3ba806937a> in <module>()
    ----> 1 myarray = linspace(0, 5, 10)
    

    NameError: name 'linspace' is not defined


The function `linspace()` is very useful. Try it changing the input parameters!

**Import style:**

You will often see code snippets that use the following lines
```Python
import numpy as np
import matplotlib.pyplot as plt
```
What's all of this import-as business? It's a way of creating a 'shortcut' to the NumPy library and the pyplot module. You will see it frequently as it is in common usage, but we prefer to keep out imports explicit. We think it helps with code readability.

**Pro tip:**

Sometimes, you'll see people importing a whole library without assigning a shortcut for it (like `from numpy import *`). This saves typing but is sloppy and can get you in trouble. Best to get into good habits from the beginning!


To learn new functions available to you, visit the [NumPy Reference](http://docs.scipy.org/doc/numpy/reference/) page. If you are a proficient `Matlab` user, there is a wiki page that should prove helpful to you: [NumPy for Matlab Users](http://wiki.scipy.org/NumPy_for_Matlab_Users)

## Variables

Python doesn't require explicitly declared variable types like C and other languages.  


```python
a = 5        #a is an integer 5
b = 'five'   #b is a string of the word 'five'
c = 5.0      #c is a floating point 5  
```


```python
type(a)
```




    int




```python
type(b)
```




    str




```python
type(c)
```




    float



Note that if you divide an integer by an integer that yields a remainder, the result will be converted to a float.  (This is *different* from the behavior in Python 2.7, beware!)

## Whitespace in Python

Python uses indents and whitespace to group statements together.  To write a short loop in C, you might use:

    for (i = 0, i < 5, i++){
       printf("Hi! \n");
    }

Python does not use curly braces like C, so the same program as above is written in Python as follows:


```python
for i in range(5):
    print("Hi \n")
```

    Hi 
    
    Hi 
    
    Hi 
    
    Hi 
    
    Hi 
    


If you have nested for-loops, there is a further indent for the inner loop.


```python
for i in range(3):
    for j in range(3):
        print(i, j)
    
    print("This statement is within the i-loop, but not the j-loop")
```

    0 0
    0 1
    0 2
    This statement is within the i-loop, but not the j-loop
    1 0
    1 1
    1 2
    This statement is within the i-loop, but not the j-loop
    2 0
    2 1
    2 2
    This statement is within the i-loop, but not the j-loop


## Slicing Arrays

In NumPy, you can look at portions of arrays in the same way as in `Matlab`, with a few extra tricks thrown in.  Let's take an array of values from 1 to 5.


```python
myvals = numpy.array([1, 2, 3, 4, 5])
myvals
```




    array([1, 2, 3, 4, 5])



Python uses a **zero-based index**, so let's look at the first and last element in the array `myvals`


```python
myvals[0], myvals[4]
```




    (1, 5)



There are 5 elements in the array `myvals`, but if we try to look at `myvals[5]`, Python will be unhappy, as `myvals[5]` is actually calling the non-existant 6th element of that array.


```python
myvals[5]
```


    ---------------------------------------------------------------------------

    IndexError                                Traceback (most recent call last)

    <ipython-input-13-6cc4d3ae83cd> in <module>()
    ----> 1 myvals[5]
    

    IndexError: index 5 is out of bounds for axis 0 with size 5


Arrays can also be 'sliced', grabbing a range of values.  Let's look at the first three elements


```python
myvals[0:3]
```




    array([1, 2, 3])



Note here, the slice is inclusive on the front end and exclusive on the back, so the above command gives us the values of `myvals[0]`, `myvals[1]` and `myvals[2]`, but not `myvals[3]`.

## Assigning Array Variables

One of the strange little quirks/features in Python that often confuses people comes up when assigning and comparing arrays of values.  Here is a quick example.  Let's start by defining a 1-D array called $a$:


```python
a = numpy.linspace(1,5,5)
```


```python
a
```




    array([ 1.,  2.,  3.,  4.,  5.])



OK, so we have an array $a$, with the values 1 through 5.  I want to make a copy of that array, called $b$, so I'll try the following:


```python
b = a
```


```python
b
```




    array([ 1.,  2.,  3.,  4.,  5.])



Great.  So $a$ has the values 1 through 5 and now so does $b$.  Now that I have a backup of $a$, I can change its values without worrying about losing data (or so I may think!).


```python
a[2] = 17
```


```python
a
```




    array([  1.,   2.,  17.,   4.,   5.])



Here, the 3rd element of $a$ has been changed to 17.  Now let's check on $b$.


```python
b
```




    array([  1.,   2.,  17.,   4.,   5.])



And that's how things go wrong!  When you use a statement like $a = b$, rather than copying all the values of $a$ into a new array called $b$, Python just creates an alias (or a pointer) called $b$ and tells it to route us to $a$.  So if we change a value in $a$ then $b$ will reflect that change (technically, this is called *assignment by reference*).  If you want to make a true copy of the array, you have to tell Python to copy every element of $a$ into a new array.  Let's call it $c$.  


```python
c = a.copy()
```

Now, we can try again to change a value in $a$ and see if the changes are also seen in $c$.  


```python
a[2] = 3
```


```python
a
```




    array([ 1.,  2.,  3.,  4.,  5.])




```python
c
```




    array([  1.,   2.,  17.,   4.,   5.])



OK, it worked!  If the difference between `a = b` and `a = b.copy()` is unclear, you should read through this again.  This issue will come back to haunt you otherwise.

## Learn More

There are a lot of resources online to learn more about using NumPy and other libraries. Just for kicks, here we use Jupyter's feature for embedding videos to point you to a short video on YouTube on using NumPy arrays.


```python
from IPython.display import YouTubeVideo
# a short video about using NumPy arrays, from Enthought
YouTubeVideo('vWkb7VahaXQ')
```





<iframe
    width="400"
    height="300"
    src="https://www.youtube.com/embed/vWkb7VahaXQ"
    frameborder="0"
    allowfullscreen
></iframe>





```python
from IPython.core.display import HTML
def css_styling():
    styles = open("../styles/custom.css", "r").read()
    return HTML(styles)
css_styling()
```




<link href='http://fonts.googleapis.com/css?family=Fenix' rel='stylesheet' type='text/css'>
<link href='http://fonts.googleapis.com/css?family=Alegreya+Sans:100,300,400,500,700,800,900,100italic,300italic,400italic,500italic,700italic,800italic,900italic' rel='stylesheet' type='text/css'>
<link href='http://fonts.googleapis.com/css?family=Source+Code+Pro:300,400' rel='stylesheet' type='text/css'>
<style>
    @font-face {
        font-family: "Computer Modern";
        src: url('http://mirrors.ctan.org/fonts/cm-unicode/fonts/otf/cmunss.otf');
    }
    div.cell{
        width:800px;
        margin-left:16% !important;
        margin-right:auto;
    }
    h1 {
        font-family: 'Alegreya Sans', sans-serif;
    }
    h2 {
        font-family: 'Fenix', serif;
    }
    h3{
		font-family: 'Fenix', serif;
        margin-top:12px;
        margin-bottom: 3px;
       }
	h4{
		font-family: 'Fenix', serif;
       }
    h5 {
        font-family: 'Alegreya Sans', sans-serif;
    }	   
    div.text_cell_render{
        font-family: 'Alegreya Sans',Computer Modern, "Helvetica Neue", Arial, Helvetica, Geneva, sans-serif;
        line-height: 135%;
        font-size: 120%;
        width:600px;
        margin-left:auto;
        margin-right:auto;
    }
    .CodeMirror{
            font-family: "Source Code Pro";
			font-size: 90%;
    }
/*    .prompt{
        display: None;
    }*/
    .text_cell_render h1 {
        font-weight: 200;
        font-size: 50pt;
		line-height: 100%;
        color:#CD2305;
        margin-bottom: 0.5em;
        margin-top: 0.5em;
        display: block;
    }	
    .text_cell_render h5 {
        font-weight: 300;
        font-size: 16pt;
        color: #CD2305;
        font-style: italic;
        margin-bottom: .5em;
        margin-top: 0.5em;
        display: block;
    }

    .warning{
        color: rgb( 240, 20, 20 )
        }  
</style>
<script>
    MathJax.Hub.Config({
                        TeX: {
                           extensions: ["AMSmath.js"]
                           },
                tex2jax: {
                    inlineMath: [ ['$','$'], ["\\(","\\)"] ],
                    displayMath: [ ['$$','$$'], ["\\[","\\]"] ]
                },
                displayAlign: 'center', // Change this to 'center' to center equations.
                "HTML-CSS": {
                    styles: {'.MathJax_Display': {"margin": 4}}
                }
        });
</script>



