---
title: D-Claw
subtitle: software for granular-fluid flows
description: D-Claw is numerical software package for simulating granular-fluid flows, such as landslides, debris flows, and lahars. D-Claw is an extension of the Clawpack codes (www.clawpack.org)
---

---

# Overview

D-Claw is numerical software for modeling granular-fluid flows.  It is built on top of Clawpack ([clawpack.org](http://www.clawpack.org)), and is an extension and generalization of the shallow-water code GeoClaw ([geoclaw.org](http://www.geoclaw.org)), which includes algorithms for general shallow earth-surface flows. D-Claw solves a more general two-phase (granular-fluid mixture) model for landslides, debris flows, and lahars. D-Claw can also be used for simulating hybrid problems that involve the interaction of granular-fluid flows with bodies of water (*eg.*, landslide-generated tsunamis, dam-breach floods, fluid or solid entrainment by inundating flows, such as debris-laden run-off or debris-entraining tsunamis). In the limit of vanishing solid concentrations, D-Claw results theoretically reduce to Geoclaw results (though details of the numerical implementation may differ). 

The documentation available with Clawpack (v5) and GeoClaw provide a general overview, however, in addition to unique features, D-Claw is built on top of older versions of Clawpack and GeoClaw (v4.x). (The older v4.x version of Clawpack is available at [github.com/clawpack/clawpack-4.x](https://github.com/clawpack/clawpack4.x). See [clawpack.org](http://www.clawpack.org) for more information.) 

Running D-Claw requires additional set-up parameters beyond those of Geoclaw. Documentation for D-Claw is currently unsatisfactory, but in progress. The basics of running it are below. See also the [geoclaw/dclaw-apps](https://github.com/geoclaw/dclaw-apps) on github for examples.


# Use

To use D-Claw as is, it is assumed that you have a unix terminal of some kind (*eg.,* linux or Mac OS...for MS Windows you are on your own, but options exist for terminal emulators). 

## Source code

The source code and latest git repository for D-Claw are available on github:

* [github.com/geoflows/D-Claw](https://github.com/geoflows/D-Claw).

A repository for applications (in progress) is also available:

* [github.com/geoflows/dclaw-apps](https://github.com/geoflows/dclaw-apps).

If you want just the source code without using git, github provides a zipped directory. Otherwise, you can clone either of the repositories with git clone as usual:

```
git clone https://github.com/geoflows/D-Claw.git
```
and
```
git clone https://github.com/geoflows/dclaw-apps.git
```
 

## Installation

"Installation" of D-Claw is essentially just setting some environment variables. Compiling and running D-Claw can then be done using make, with a Makefile in an application directory.  

#### environment variables

The environment variable $CLAW should be set to point to the parent directory of the D-Claw source code. For bash shells:
```
export CLAW=/path/to/D-Claw
```

TIP: if you are using multiple versions of Clawpack (*eg.,* Clawpack 5.x or GeoClaw and D-Claw), it is advisable to stay mindful of how $CLAW is set, particularly if you are developing/testing code. It is easy to inadvertently compile the wrong code, or fail to incorporate code-changes if $CLAW is set incorrectly. Packages such as the [environment modules](http://modules.sourceforge.net/) package, which can dynamically set or change your environment under a given shell to ensure that you have a compatible set of paths/versions of software (*eg.*, $PATH, $CLAW, $PYTHONPATH, $MATLABPATH), are very useful.

(Alternatively, one could modify Makefiles if you want to use an environment variable other than $CLAW, such as $DCLAW. However, this can add complications to your path hierarchy for python and matlab when switching between different $CLAWs, as there are modified routines with the same file name...and so it is not recommended without something like environmental modules.)   

#### python
In order to use the python set-up scripts that come with D-Claw, you should include $CLAW/python/ in your $PYTHONPATH:
```
export PYTHONPATH=$CLAW/python:$PYTHONPATH
```
You can also then import some auxiliary tools from the D-Claw python folder if you find them useful:  
```
python> import dclaw
python> import dclaw.topotools as dt 
```
See the comments in these modules for more information.

NOTE: The current version of D-Claw is compatible with python v2.x. If you prefer python version v3.x, see below about developing and contributing to D-Claw.    

#### matlab
If you are going to use matlab plotting with D-Claw, your $MATLABPATH should include $CLAW/matlabgeo, above but in addition to $CLAW/matlab:
```
export MATLABPATH=$CLAW/matlabgeo:$CLAW/matlab:$MATLABPATH
```

## Running a D-Claw simulation

Running D-Claw applications is very similar to running applications in Clawpack v5.x. See the documentation at [clawpack.org](http://www.clawpack.org). See also example applications that include necessary files in the application repository, [github/geoflows/dlcaw-apps](https://github.com/geoflows/dclaw-apps).

The D-Claw executable is made entirely from Fortran source code. However, python is used to set-up and initialize a given D-Claw run, and optionally plot the results. Therefore, ordinary use of D-Claw requires interaction with python scripts and the make program. Modifying Fortran source code is only necessary if you are developing new features or debugging.

#### directory and file structure

In summary, a working application directory (it is recommended that this reside away from your D-Claw source code) for a single simulation should contain:
* a Makefile
* python initialization scripts (setrun.py)
* optional plotting routines (python or matlab)
* optional pre-processing routines for topography or other data requirements

Setting up a given simulation essentially amounts to modifying the routine setrun.py, to set runtime parameters, initial conditions, and required input data (*eg.,* topography DEMs). Plotting can be done with python or matlab scripts, which are included with D-Claw and application examples (*see* *eg.,* [geoclaw/dclaw-apps](https://github.com/geoclaw/dclaw-apps) on github). See more info on plotting below, or in example applications.

#### make 


```
pwd

/path/to/myapplication/
```

To recompile all of the Fortran source code from scratch and create a new executable ("xgeoclaw"), issue:

```
make new
```
 
The following make commands will take into account dependency changes to determine the sequence of prerequisite steps that need to be taken. To make or retain an up-to-date executable, issue: 
```
make .exe
```

To run the executable and produce output (the output files will be placed in a subdirectory indicated in the Makefile), issue:
```
make .output
```

To produce plots using python and the setplot.py routine in your application directory, issue: 
```
make .plots
```

Note that each one of the above steps depends on the previous steps if source code or parameters have changed. So, for instance, "make .plots" will recompile source code, rerun the executable to produce new output, and finally produce new plots if the source code has changed. If only setrun.py has been modified, it will re-run the existing executable to produce new output. If nothing has changed, make will indicate that nothing needs to be done.

**WARNING: if the commands `make .output` or `make .plots` create a new run, previously made output files in the `_output` directory will be deleted. If you have changed source code or runtime parameters, but wish to keep your old output for comparison or debugging, take necessary action to save your old files. The Makefile can also be modified to specify your output directory name.** 

To produce plots with python without checking dependencies, without the risk of deleting wanted output files, you can issue:
```
make plots
```
 (without the ".").

## Plotting results
#### matlab

Matlab can be used to plot D-Claw output. From the output directory, use
```
matlab> plotclaw2
```
then follow the interactive menu to produce plots for each frame.

TIP: For a given application, it is useful to relocate some of the m-files (*eg.,* afterframe.m, setplot2.m, setprob.m, beforeframe.m etc.) included with D-Claw to your local working application directory, where you can modify them to suit you present purposes without modifying your D-Claw source files (note that files in $CLAW/matlabgeo should take precedent over files of the same name in $CLAW/matlab. Unifying these directories is planned in the future, but they currently coexist so that the D-Claw source code can be used to run Geoclaw v4.x applications...more about that in the future). 

TIP: Note that D-Claw output and these locally modified .m-files (in your application directory) must both be located by matlab, but they are not usually in the same directory. For instance, if the output sub-directory is the current directory where matlab is running, *ie.,*
```
matlab> pwd
/path/to/myapplication/_output
```
then you can issue,
```
matlab> addpath ../
```
to get the local .m-files in the output's parent directory, myapplication/, to the top of your path (*ie.* Matlab will add the absolute path for ../ to the top of your path for the current session).

NOTE: you could alternatively place your local m-files in the output directory...but this is not recommended if you want your local m-files to be part of a repository, as the output directory is best ignored by git, as it is with the applications in the [dclaw-apps](https://github.com/geoflows/dclaw-apps).

More information on plotting can be found in the [dclaw-apps](https://github.com/geoflows/dclaw-apps) repository.


#### python

Python can alternatively be used to produce mapview 2D or 1D cross-sectional plots, as describe above, with:
```
make .plots
```
or
```
make plots
```
Make uses setplot.py and matplotlib.  Modify the routine setplot.py to your needs. 

See [clawpack.org](http://www.clawpack.org), [github/clawpack/visclaw](https://gihub.com/clawpack/visclaw) and [github/geoflows/dclaw-apps](https://github.com/geoflows/dclaw-apps) for more information.


## Development

* If you would like to make contributions to D-Claw or dclaw-apps, please follow the development workflow used for Clawpack, described at [www.clawpack.org/developers](http://www.clawpack.org/developers.html). In summary, please fork the repositories to your own github account, and issue pull requests on a feature branch to github/geoflows, *eg,*:

```
git clone https://github.com/geoflows/D-Claw.git
cd D-Claw
git remote add username https://github.com/username/D-Claw.git
```
or if you have ssh keys and want to avoid typing your password when you push to github:

```
git remote add username git@github.com:username/D-Claw.git
```
These settings can be modified in your local working repository at anytime with `git remote set-url`.

* Develop in a branch other than master:
```
git checkout -b my_branch
```
And then push to your repository on github:
```
git push username my_branch
```
* Issue pull requests from your branch and repository on github.com (username/D-Claw) to contribute features or fixes to the D-Claw master branch at geoflows/D-Claw. 

* Update your master branches from geoflows/D-Claw:
```
git pull origin master
```
and then 
```
git push username master
```
to update your git remote. It is recommended that you do not commit to your own master branches, so that your master branches are easily updated from the geoflows repository.

If you prefer, rename origin to something easy to remember ("geoflows" or "upstream" or similar):
```
git remote rename origin geoflows
```

## License

D-Claw inherits the Clawpack licenses and user agreeements. 