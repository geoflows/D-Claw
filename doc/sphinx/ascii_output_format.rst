
.. _ascii_output_format:

Format of ASCII output files
----------------------------

At each output time, two files are created:

* fort.tNNNN  small file with some metadata for frame N

* fort.qNNNN  file with the solution on all grids for frame N


`fort.t` file
-------------

Sample file in 2 dimensions::

        0.00000000E+00    time
        4                 meqn
        2                 ngrids
        3                 naux
        2                 ndim

`fort.q` file
-------------

This file contains all solution value from each of the `ngrids` grids.
Each set of grid data is preceeded by a header of the form::

    1                 grid_number
    1                 AMR_level
  200                 mx
  200                 my
    0.00000000E+00    xlow
    0.00000000E+00    ylow
    0.87500000E-01    dx
    0.87500000E-01    dy

followed by `mx * my` lines, each having `meqn` values of q, in the `(i,j)`
grid cell, created by this loop in `out2.f`::

          do j=1,my
            do i=1,mx
              do m=1,meqn
    c            # exponents with more than 2 digits cause problems 
    c            # reset tiny values to zero:
                 if (dabs(q(i,j,m)) .lt. 1d-99) q(i,j,m) = 0.d0
                 enddo
              write(50,1005) (q(i,j,m), m=1,meqn)
     1005     format(4e26.16)

This grid is followed immediately by the header for the next grid.

