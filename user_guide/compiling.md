# D-Claw Compiling tips and tricks

Katy Barnhart


## issues with call.i

  On linux two files called call.i are symlinks. If you clone on a mac or a pc, you will probably get an error that looks like:

  ```bash
  Error: Non-numeric character in statement label at (1)
  call.i:1:2:

  Error: Invalid character in name at (1)
  ```

  To address this, rebuild the symlinks

  From the D-Claw directory

  ```bash
  cd geoclaw/2d/lib
  rm call.i
  ln -s ../../../amrclaw/2d/lib/call.i call.i
  cd ../lib_dig
  rm call.i
  ln -s ../../../amrclaw/2d/lib/call.i call.i
  ```


## issues with type missmatch causing compile errors

  If you have fortran compiler > v8 (I think) you will probably get the error:

  ```bash
  Error: Type mismatch in argument ‘iflags’ at (1); passed REAL(8) to INTEGER(1)
  ```

  referring to spest.f in amrclaw.

  I do not get this with GNU Fortran (GCC) 8.3.0 20190222 but get it with GNU Fortran (GCC) 11.2.0

  add `-fallow-argument-mismatch` to $FFLAGS
