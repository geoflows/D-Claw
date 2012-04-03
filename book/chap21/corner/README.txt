

begin_html [use: jsMath] [use: doc/doc.css]
<!--   For a more readable version of this file, execute
                  unix>  make htmls
       in this directory and then point your browser to README.html
     --------------------------------------------------------------  -->

<h2>
CLAWPACK Sample Code
</h2>

2-d acoustics with variable coefficients.

Two materials with a single interface, specified in [code: fdisc.f]
(which is called from [code: cellave.f] to compute the fraction of
each grid cell lying in the left and right states).

<img src="interface.png" width=300> Figure created by [code: makegridfig.py]

The sound speed and impedance are stored in the aux array, specified in
[code: setaux.f].
For cells cut by the interface, the arithmetic average of rho and harmonic
average of the bulk modulus are used to determine these parameters for this
cell.
    
Example [book/chap21/corner]
to accompany the book <br> &nbsp;&nbsp;
  [www.clawpack.org/book.html Finite Volume Methods for Hyperbolic Problems]
  by R. J. LeVeque.

Converted to [www.clawpack.org Clawpack 4.5] form in 2011.
        
<h4>
Plots of results
</h4>
After running this code and creating plots via "make .plots", you should be
able to view the plots in [link: _plots/_PlotIndex.html].


<h4>
Fortran files
</h4>

<dl>
<dt>[code: Makefile]
<dd> Determines which version of fortran files
are used when compiling the code with make and specifies where output and
plots should be directed.  Type "make .help" at the Unix prompt for options.

<dt>[code: driver.f]
<dd>
The driver routine allocates storage and then calls the main Clawpack
routine.

<dt>[code: setprob.f]
<dd>
A routine by this name is called by the library routine
[clawcode: clawpack/2d/lib/claw2ez.f]
and is generally used to set any values needed for the specific problem
being solved.
    
    
<dt>[code: rpn2acv.f]
<dd> Normal Riemann solver (normal to cell interface).

<dt>[code: rpt2acv.f]
<dd> Transverse Riemann solver.

         
<dt>[code: qinit.f]
<dd>
This subroutine sets the initial data at time $t=0$.

<dt>[code: setaux.f]
<dd>
Sets the aux arrays.

<dt>[code: fdisc.f]
<dd>
Determines the interface in material properties.


</dl>

<h4>
Python files
</h4>
<dl>

<dt>[code: setrun.py]
<dd> This file contains a function that
specifies what run-time parameters will be used.

Some parameters that you might want to modify are described in the
[www.clawpack.org/doc.html documentation].

<dt>[code: setplot.py]
<dd> This file contains a function that
specifies what plots will be done and
sets various plotting parameters.

</dl>


<h4>
Data files
</h4>
<font color='red'>Warning:</font> These files are generally changed
when setting up a run, usually in [code: setrun.py].

<dl>
<dt>[code: claw.data]
<dd> This file contains a number of
parameter values that are used by CLAWPACK.
The values in this file are read by the library routine
[clawcode: clawpack/2d/lib/claw2ez.f].


<dt> [code: setprob.data]
<dd> This file may contain various
parameters used in setting the initial conditions or otherwise setting up
the problem.


</dl>

<h4>Library routines</h3>

In addition to the Fortran routines in this library, several library
routines from [claw:clawpack/2d/lib] are used.  See the [code: Makefile]
to determine which ones are used.

end_html

    
