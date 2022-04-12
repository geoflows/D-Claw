# D-Claw Compiling tips and tricks

Katy Barnhart


## Compilation fails: call.i
  In the repo, two call.i files are intended as symlinks. They show up as text documents on Mac or Windows and if they are not fixed, you will probably get the following error message.

  ```bash
  Error: Non-numeric character in statement label at (1)
  call.i:1:2:

  Error: Invalid character in name at (1)
  ```

  If on Mac or Windows, need to remove these files and rebuild the symlinks.

  ```
  cd D-Claw/geoclaw/2d/lib
  rm call.i
  ln -s ../../../amrclaw/2d/lib/call.i call.i
  cd ../lib_dig
  rm call.i
  ln -s ../../../amrclaw/2d/lib/call.i call.i
  ```

## Compilation fails: type mismatch error 

  If you have fortran compiler > v8 (I think) you will probably get the error:

  ```bash
  Error: Type mismatch in argument ‘iflags’ at (1); passed REAL(8) to INTEGER(1)
  ```

  referring to spest.f in amrclaw.

  I do not get this with GNU Fortran (GCC) 8.3.0 20190222 but get it with GNU Fortran (GCC) 11.2.0

  add `-fallow-argument-mismatch` to $FFLAGS
