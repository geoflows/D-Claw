# DCLAW Windows installation instructions and notes
Katy Barnhart

System description: Windows 10 machine.

Additional goals: Be able to edit files on Windows side and have them run on the linux side. Be able to inspect and post-process results natively on the windows side.

1. Install WSL
2. Install UBUNTU 18.04

On Sept 6, 2022 we had a bunch of permissions issues related to wsl-windows file permissions.

https://devblogs.microsoft.com/commandline/chmod-chown-wsl-improvements/

Second answer here worked:
https://askubuntu.com/questions/1115564/wsl-ubuntu-distro-how-to-solve-operation-not-permitted-on-cloning-repository Quoting from this answer:
On WSL edit `/etc/wsl.conf` (create it if it doesn't exist). Add the following:

```bash
[automount]
options = "metadata"
```
  Then either:

  Reboot Windows
  Exit any WSL sessions, run `wsl --shutdown` from PowerShell or CMD, and start WSL again
  Exit your only session, terminate it with `wsl --terminate <distroName>`, and start it again,

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

  c. Then install miniconda and ensure conda is updated. Might need to remove the slash (seems to depend based on platform).
  ```bash
  bash /Miniconda3-latest-Linux-x86_64.sh
  conda update -n base -c defaults conda
  ```

5. Create compute environment (this is probably where to start if you are on a mac)

  a. Create python environment with necessary dependencies.
  ```bash
  cd D-Claw
  conda env create --file=environment.yml
  conda activate dclaw
  ```

  b. Compile the python parts of D-Claw (this will help access the D-Claw python code).

  ```bash
  cd python
  pip install -e .
  ```

6. Within the D-Claw repository configure paths correctly. One should expect to compile D-Claw for each simulation evaluation. Thus "installation" takes the form of setting environment variables and paths correctly.

  a. Inspect `setenv.py` and modify if necessary. I found no changes necessary.

  b. Run the following to set environments correctly.

  You'll need to go up one directory level.

  ```bash
  python setenv.py
  source setenv.bash
  ```
  You'll get a message like:

  ------------------------------------------------------------
Full path to claw directory should be:
      $CLAW =  /Users/krbarnhart/krbarnhart/source/D-Claw
------------------------------------------------------------
The files setenv.csh and setenv.bash contain the appropriate
commands to set environment variables for csh or bash shells
  and also some aliases you may find convenient
------------------------------------------------------------

This means TODO.

  EVERY TIME YOU RUN DCLAW, THESE ENVIRONMENT VARIABLES MUST BE SET. THIS MEANS YOU EITHER NEED TO COPY THIS TEXT INTO AN ENVIRONMENT FILE (E.G., FOR SLURM) OR YOU NEED TO NAVIGATE HERE AND SOURCE SETENV.BASH
TODO> ALSO NEED TO activeate the conda env.

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

7. Test based on running an example. I used the USGS gate release example. Following the instructions in  ``dclaw-apps/USGSFlume/gate_release_example\readme.md`` a minimal approach to this is:

  ```bash
  python setinit.py
  make new
  make .plots
  ```

  This worked for me. It took maybe 10 minutes. If you have compiler errors, look at compiling.md.
  Sept 2022: will probably get a compiler error related to type mismatch.

  Presuming you have a latex distribution installed you can also compile a PDF of the figures with

  ```bash
  cd _plots
  pdflatex plots.tex
  ```

## Other notes

Visual Studio Code and MobaXTerm seem like the best choices if you want to "remote" into your WSL installation. MobaXTerm can also serve as your display viewer if you, say, want to run paraview on WSL.
