# DCLAW Windows installation instructions and notes
Katy Barnhart

System description: Windows 10 machine.

Additional goals: Be able to edit files on Windows side and have them run on the linux side. Be able to inspect and post-process results natively on the windows side.

1. Install WSL
2. Install UBUNTU 18.04
3. In a WSL terminal, create a symlink to the folders _on the Windows side_ where you want to have your source files. This will allow you to edit with native Windows apps and run on WSL.

  For example I made a directory called `source`.

  ```bash
  ln -s /mnt/c/Users/krbarnhart/data/source source
  ```

  Critical Note. From the WSL terminal, into the symlinked folder, download the following source distributions using git.

  a. D-CLAW https://github.com/geoflows/D-Claw

  b. D-CLAW apps https://github.com/geoflows/dclaw-apps

  ```bash
  git clone https://github.com/geoflows/D-Claw.git
  git clone https://github.com/geoflows/dclaw-apps.git
  ```

  The only way I found that worked was to clone DIRECTLY to the WSL file structure. If you clone the source to the windows side you may get fortran compiler errors. If you clone to WSL without setting up the symlinks then you can't edit directly on windows (as easily).

  One thought I had was that this might be a fortran compiler fighting with windows line endings issue (CRLF vs LF). So I tried cloning to the Windows file structure with git configured to use either default local line endings (CRLF) OR source line endings (presumably LF), and under these circumstances got the same compiler error at step (X) from below:
  ```
  gfortran -c -I/mnt/c/Users/krbarnhart/source/D-Claw/amrclaw/2d/lib /mnt/c/Users/krbarnhart/source/D-Claw/geoclaw/2d/lib_dig/setaux_geo.f -o /mnt/c/Users/krbarnhart/source/D-Claw/geoclaw/2d/lib_dig/setaux_geo.o
  call.i:1:1:

  Error: Non-numeric character in statement label at (1)
  call.i:1:2:

  Error: Invalid character in name at (1)
  /mnt/c/Users/krbarnhart/source/D-Claw/util/Makefile.common:71: recipe for target '/mnt/c/Users/krbarnhart/source/D-Claw/geoclaw/2d/lib_dig/setaux_geo.o' failed
  make[1]: *** [/mnt/c/Users/krbarnhart/source/D-Claw/geoclaw/2d/lib_dig/setaux_geo.o] Error 1
  make[1]: Leaving directory '/mnt/c/Users/krbarnhart/source/dclaw-apps/USGSFlume/gate_release_example'
  ```

4. On WSL do additional steps of installation.

  a. Added gfortran and make to WSL
  ```bash
  sudo apt install make
  sudo apt install gfortran
  ```

  b. Add miniconda to WSL. First download from anaconda.com with wget.
  ```bash
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86.sh
  ```

  c. Then install miniconda and ensure conda is updated.
  ```bash
  bash /Miniconda3-latest-Linux-x86.sh
  conda update -n base -c defaults conda
  ```

  d. Create python 2.7 environment.
  ```bash
  conda create -n py27 python=2.7
  ```

  f. Add core python dependencies and clean up.

  ```bash
  conda activate py27
  conda install numpy scipy matplotlib
  conda clean --all
  ```

  scipy seems to be a requirement but was not listed as such.

  There are other potentially relevant packages to install (based on casting about the python code, but running examples has not yet necessitated these). Specifically:
    - petsc4py
    - mpi4py
    - h5py
    - ipython

  As I gain an understanding about whether/if these are needed, I'll update this to describe for what.

  I also installed ipython and jupyter because I find them useful.

  ```bash
  conda activate py27
  conda install ipython jupyter
  conda clean --all
  ```
5. Within the D-Claw repository configure paths correctly. One should expect to compile D-Claw for each simulation evaluation. Thus "installation" takes the form of setting environment variables and paths correctly.

  a. Inspect `setenv.py` and modify if necessary. I found no changes necessary.

  b. Run the following to set environments correctly.

  ```bash
  python setenv.py
  source setenv.bash
  ```

  On my WSL install the contents of `setenv.bash` is as follows:
  ```
  export CLAW='/mnt/c/Users/krbarnhart/source/D-Claw'
  export FC='gfortran'

  export MATLABPATH='/mnt/c/Users/krbarnhart/source/D-Claw/matlab'

  export PYTHONPATH="/mnt/c/Users/krbarnhart/source/D-Claw/python:${PYTHONPATH}"
  export IPYTHONDIR='/mnt/c/Users/krbarnhart/source/D-Claw/python/ipythondir'
  if [ -z "${LD_LIBRARY_PATH}" ]; then
      export LD_LIBRARY_PATH="/mnt/c/Users/krbarnhart/source/D-Claw/lib"
  else
      export LD_LIBRARY_PATH="/mnt/c/Users/krbarnhart/source/D-Claw/lib:${LD_LIBRARY_PATH}"
  fi
  alias ipyclaw='ipython -profile claw'
  alias clawserver='xterm -e python $CLAW/python/startserver.py &'
```

6. Test based on running an example. I used the USGS gate release example. Following the instructions in  ``dclaw-apps/USGSFlume/gate_release_example\readme.md`` a minimal approach to this is:

  ```bash
  python setinit.py
  make .plots
  ```

  This worked for me. It took maybe 10 minutes.

  Presuming you have a latex distribution installed you can also compile a PDF of the figures with

  ```bash
  cd _plots
  pdflatex plots.tex
  ```

## Other notes

Visual Studio Code and MobaXTerm seem like the best choices if you want to "remote" into your WSL installation. MobaXTerm can also serve as your display viewer if you, say, want to run paraview on WSL.
