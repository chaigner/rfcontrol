rfcontrol [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.14808.svg)](http://dx.doi.org/10.5281/zenodo.31187)
========

This repository contains a MATLAB implementation of the approach for computing RF pulses in magnetic resonance imaging described in [1], based on a trust-region CG-Newton method for optimal control of the full time-dependent Bloch equation. Test scripts are provided to reproduce the numerical optimization and simulation results in the paper [1].

##### Authors:
- Christoph Aigner    (<christoph.aigner@tugraz.at>)
- Christian Clason    (<christian.clason@uni-due.de>)
- Armin Rund          (<armin.rund@uni-graz.at>)

Contents
--------

##### Test scripts (run these):
    test_single.m          test script for single-slice example
    test_multi.m           test script for multi-slice example

##### Routines called by the test scripts:
    tr_newton.m            implements trust-region Newton method
    tr_cg.m                implements trust-region conjugate gradient iteration
    objfun.m               computes functional value, gradient
    applyHess.m            computes action of Hessian
    cn_bloch.m             solves Bloch equation using Crank-Nicolson scheme
    cn_adjoint.m           solves adjoint equation using adjoint CN scheme
    plot_results.m         plots optimized pulse, magnetization
    
##### Data files used by the test scripts:
    z_grad_thk2_dt5.mat    binary file containing GRE gradient 


Dependencies
------------

These routines were tested under MATLAB R2013b -- R2014b under Linux (x86_64), OS X (10.10) and Windows, but should also run under older versions. The `Parallel Toolbox` is used to accelerate the solution of the differential equations via `parfor` loops; for MATLAB versions prior to R2013b, make sure to start a `matlabpool` before running the scripts. If the toolbox is not installed, these loops should automatically be executed in serial. For older versions of MATLAB not implementing the `parfor` keyword (pre 2008a), simply replace all occurences of `parfor` by `for` in `cn_bloch.m`, `cn_adjoint.m` and `applyHess.m`.

License
-------

This software is published under GNU GPLv3. 
In particular, all source code is provided "as is" without warranty of any kind, either expressed or implied. 
For details, see the attached LICENSE.

Reference
---------

[1] C. S. Aigner, C. Clason, A. Rund and R. Stollberger: <br/>
&nbsp;&nbsp;&nbsp;&nbsp;*Efficient high-resolution RF pulse design applied to simultaneous 
multi-slice excitation*, <br/>
&nbsp;&nbsp;&nbsp;&nbsp;[Journal of Magnetic Resonance, in press, 2015](http://dx.doi.org/10.1016/j.jmr.2015.11.013)
([Preprint](http://math.uni-graz.at/mobis/publications/SFB-Report-2015-001.pdf)).<br/>

##### bibtex:

Please cite this work as 

    @article{ rfcontrol,
       author      = {Aigner, Christoph Stefan and Clason, Christian and Rund, Armin and Stollberger, Rudolf},
       title       = {Efficient high-resolution {RF} pulse design applied to simultaneous multi-slice excitation},
       journal     = {Journal of Magnetic Resonance},
       volume      = {263},
       pages       = {33--44},
       year        = {2016},
       doi         = {10.1016/j.jmr.2015.11.013}
    }

If you wish to cite this implementation specifically, you can do this as

    @misc{ rfcontrol,
      author = {Aigner, Christoph Stefan and Clason, Christian and Rund, Armin},
      title  = {rfcontrol},
      year   = {2015},
      doi    = {10.5281/zenodo.31187}
    }
