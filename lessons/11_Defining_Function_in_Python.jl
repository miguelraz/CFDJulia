Text provided under a Creative Commons Attribution license, CC-BY.  All code is made available under the FSF-approved BSD-3 license.  (c) Lorena A. Barba, Gilbert F. Forsyth 2017. Thanks to NSF for support via CAREER award #1149784.
[@LorenaABarba](https://twitter.com/LorenaABarba)

12 steps to Navier–Stokes
=====
***

This lesson complements the first interactive module of the online [CFD Python](https://github.com/barbagroup/CFDPython) class, by Prof. Lorena A. Barba, called **12 Steps to Navier–Stokes.** The interactive module starts with simple exercises in 1D that at first use little of the power of Python. We now present some new ways of doing the same things that are more efficient and produce prettier code.

This lesson was written with BU graduate student Gilbert Forsyth.


Defining Functions in Python 
----

In steps 1 through 8, we wrote Python code that is meant to run from top to bottom.  We were able to reuse code (to great effect!) by copying and pasting, to incrementally build a solver for the Burgers' equation. But moving forward there are more efficient ways to write our Python codes.  In this lesson, we are going to introduce *function definitions*, which will allow us more flexibility in reusing and also in organizing our code.  

We'll begin with a trivial example: a function which adds two numbers.  

To create a function in Python, we start with the following:

    def simpleadd(a,b):

This statement creates a function called `simpleadd` which takes two inputs, `a` and `b`. Let's execute this definition code.


```python
def simpleadd(a, b):
    return a+b
```

The `return` statement tells Python what data to return in response to being called.  Now we can try calling our `simpleadd` function:


```python
simpleadd(3, 4)
```




    7



Of course, there can be much more happening between the `def` line and the `return` line.  In this way, one can build code in a *modular* way. Let's try a function which returns the `n`-th number in the Fibonacci sequence.  


```python
def fibonacci(n):
    a, b = 0, 1
    for i in range(n):
        a, b = b, a + b
    return a
    
```


```python
fibonacci(7)
```




    13



Once defined, the function `fibonacci` can be called like any of the built-in Python functions that we've already used. For exmaple, we might want to print out the Fibonacci sequence up through the `n`-th value:


```python
for n in range(10):
    print(fibonacci(n))
```

    0
    1
    1
    2
    3
    5
    8
    13
    21
    34


We will use the capacity of defining our own functions in Python to help us build code that is easier to reuse, easier to maintain, easier to share!

##### Exercise

(Pending.)

Learn more
-----
***

Remember our short detour on using [array operations with NumPy](./07_Step_5.ipynb)?

Well, there are a few more ways to make your scientific codes in Python run faster. We recommend the article on the Technical Discovery blog about [Speeding Up Python](http://technicaldiscovery.blogspot.com/2011/06/speeding-up-python-numpy-cython-and.html) (June 20, 2011), which talks about NumPy, Cython and Weave. It uses as example the Laplace equation (which we will solve in [Step 9](./12_Step_9.ipynb)) and makes neat use of defined functions.

But a recent new way to get fast Python codes is [Numba](http://numba.pydata.org). We'll learn a bit about that after we finish the **12 steps to Navier–Stokes**.

There are many exciting things happening in the world of high-performance Python right now!

***


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




> (The cell above executes the style for this notebook.)
