---
title: D-Claw
subtitle: software for granular-fluid flows
description: D-Claw is numerical software package for simulating granular-fluid flows, such as landslides, debris flows, and lahars. D-Claw is an extension of the Clawpack codes (www.clawpack.org)
---

---

# Summary

D-Claw is numerical software for modeling granular-fluid flows.  It is built on top of Clawpack ([clawpack.org](http://www.clawpack.org)), and is an extension and generalization of the shallow-water code GeoClaw ([geoclaw.org](http://www.geoclaw.org)), which includes algorithms for general shallow earth-surface flows. D-Claw solves a more general two-phase (granular-fluid mixture) model for landslides, debris flows, and lahars. D-Claw can also be used for simulating hybrid problems that involve the interaction of granular-fluid flows with bodies of water (*eg.*, landslide-generated tsunamis, dam-breach floods, fluid or solid entrainment by inundating flows, such as debris-laden run-off or debris-entraining tsunamis). In the limit of vanishing solid concentrations, D-Claw theoretically reduces to Geoclaw. 

The documentation available with Clawpack (v5) and GeoClaw provide a general overview, however, in addition to unique features, D-Claw is built on top of older versions of Clawpack and GeoClaw (v4.x). (The older v4.x version of Clawpack is available at [github.com/clawpack/clawpack-4.x](https://github.com/clawpack/clawpack4.x). See [clawpack.org](http://www.clawpack.org) for more information.) 

Running D-Claw requires additional set-up parameters. Documentation for D-Claw is currently unsatisfactory, but in progress. Some tips on running it are below. See also the [geoclaw/dclaw-apps](https://github.com/geoclaw/dclaw-apps) on github for examples.


# Documentation (in progress)

To use D-Claw as is, it is assumed that you have a unix terminal of some kind (*eg.,* linux or Mac OS...for MS Windows you are on your own). 

## Source code

The source code and latest git repository for D-Claw are available on github:

* [github.com/geoflows/D-Claw](https://github.com/geoflows/D-Claw).

A repository for applications (in progress) is also available:

* [github.com/geoflows/dclaw-apps](https://github.com/geoflows/dclaw-apps).
 

## Installation

"Installation" of D-Claw is essentially just setting some environment variables. Compiling and running D-Claw can then be done using make, with a Makefile in an application directory.  

#### environment variables

The environment variable $CLAW should be set to point to the parent directory of the D-Claw source code. For bash shells:
```
export CLAW=/path/to/D-Claw
```
(Alternatively, one could modify Makefiles if you want to use an environment variable other than $CLAW, such as $DCLAW.)  

TIP: if you are using multiple versions of Clawpack (*eg.,* Clawpack 5.x or GeoClaw and D-Claw), it is advisable to stay mindful of how $CLAW is set, particularly if you are developing/testing code. It is easy to inadvertently compile the wrong code, or fail to incorporate code-changes if $CLAW is set incorrectly. Packages such as the [environment modules](http://modules.sourceforge.net/) package, which can dynamically set or change your environment under a given shell to ensure that you have a compatible set of paths/versions of software (*eg.*, $PATH, $CLAW, $PYTHONPATH, $MATLABPATH), are very useful. 

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
NOTE: The current version of D-Claw is compatible with python v2.x. If you prefer python version v3.x, see below about developing and contributing to D-Claw.    

#### matlab
If you are going to use matlab plotting with D-Claw, your $MATLABPATH should include $CLAW/matlabgeo, above but in addition to $CLAW/matlab:
```
export MATLABPATH=$CLAW/matlabgeo:$CLAW/matlab:$MATLABPATH
```

## Running a D-Claw simulation

Running D-Claw applications is very similar to running applications in Clawpack v5.x. See the documentation at [clawpack.org](http://www.clawpack.org). 

#### file structure

In summary, a working application directory (it is recommended that this reside away from your D-Claw source code) for a single simulation should contain:
* a Makefile
* python initialization scripts (setrun.py)
* optional plotting routines (python or matlab)
* optional pre-processing routines for topography or other data requirements

Setting up a given simulation essentially amounts to modifying the routine setrun.py, to set runtime parameters, initial conditions, and required input data (*eg.,* topography DEMs). Plotting can be done with python or matlab scripts, which are included with D-Claw and application examples (*see* *eg.,* [geoclaw/dclaw-apps](https://github.com/geoclaw/dclaw-apps) on github) See more below.

Unfortunately documentation is still lacking for the unique attributes of D-Claw. However, examples that can be modified can be found in the [geoflows/dclaw-apps](https://github.com/geoflows/dclaw-apps) repository.

#### make 

The program make is used to compile, run, and optionally plot D-Claw output, all from your application directory:

```
pwd

/path/to/myapplication/
```

To recompile all source code from scratch and create a new executable ("xgeoclaw"), issue:

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

Note that each one of the above steps depends on the previous steps if source code or parameters have changed. So, for instance, "make .plots" will recompile source code, rerun the executable to produce new output, and finally produce new plots if the source code has changed. If nothing has changed, make will indicate that nothing needs to be done.

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
to get the local .m-files in the output's parent directory, myapplication/, to the top of your path (*ie.* Matlab will add the absolute path for ../ to the top of your path).

#### python

Python can alternatively be used to produce mapview 2d plots, using setplot.py and matplotlib. See [clawpack.org](http://www.clawpack.org) and [github/clawpack/visclaw](https://gihub.com/clawpack/visclaw) for more information about plotting with python. Note that Clawpack's v5.x python libraries may not be compatible with D-Claw v4.x output. 

## Development

If you would like to make contributions to D-Claw or dclaw-apps, please follow the development workflow used for Clawpack, described at [www.clawpack.org/developers](http://www.clawpack.org/developers.html). In summary, please fork the repositories to your own github account, and issue pull requests on a feature branch to github/geoflows, *eg,*:

```
git clone git://github.com/geoflows/D-Claw.git
cd D-Claw
git remote add username htpps://github.com/username/D-Claw.git
```
or if you have ssh keys and want to avoid typing your password when you push to github:

```
git remote add username git@github.com:username/D-Claw.git
```
Develop in a branch other than master:
```
git checkout -b my_branch
```
And then push to your repository:
```
git push username my_branch
```
Issue pull requests to geoflows/D-Claw from your repository to contribute to D-Claw. Update your master branches from geoflows/D-Claw:
```
git pull origin master
```
If you prefer, rename origin to something easy to remember ("geoflows" or "upstream" or similar):
```
git remote rename origin geoflows
```

## License

D-Claw inherits the Clawpack licenses and user agreeements. 