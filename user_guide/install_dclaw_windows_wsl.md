# DCLAW Windows installation instructions and notes
Katy Barnhart

System description: Windows 10 machine.

Additional goals: Be able to edit files on Windows side and have them run on the linux side. Be able to inspect and post-process results natively on the windows side.

1. Install WSL
2. Install UBUNTU 18.04
3. In a WSL terminal, create a symlink to the folders _on the Windows side_ where you want to have your source files. This will allow you to edit with native Windows apps and run on WSL.

  For example I made a directory called `source`.

  First make this on the windows side (wherever you'd like it located) and then create a simlink for ease of access while in WSL.

  ```bash
  ln -s /mnt/c/Users/username/source source
  ```

  Critical Note. From the WSL terminal, into the symlinked folder, download the following source distributions using git.

  a. D-CLAW https://github.com/geoflows/D-Claw

  b. D-CLAW apps https://github.com/geoflows/dclaw-apps

  ```bash
  git clone https://github.com/geoflows/D-Claw.git
  git clone https://github.com/geoflows/dclaw-apps.git
  ```

  The only way I found that worked was to clone DIRECTLY to the WSL file structure. If you clone the source to the windows side you may get fortran compiler errors. If you clone to WSL without setting up the symlinks then you can't edit directly on windows (as easily).

  See page on compiling for issues there.

4. On WSL do additional steps of installation.

  a. Added gfortran and make to WSL
  ```bash
  sudo apt install make
  sudo apt install gfortran
  ```

  b. Add miniconda to WSL. First download from anaconda.com with wget.
  ```bash
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  ```

  c. Then install miniconda and ensure conda is updated.
  ```bash
  bash /Miniconda3-latest-Linux-x86_64.sh
  conda update -n base -c defaults conda
  ```

  d. Create python environment with necessary dependencies.
  ```bash
  cd D-Claw
  conda env create --file=environment.yml
  conda activate dclaw
  ```

  e. Compile the python parts of D-Claw (this will help access the D-Claw python code).

  ```bash
  cd D-Claw\python.
  pip install -e .
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
