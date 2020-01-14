---
title: D-Claw
subtitle: software for granular-fluid flows
description: D-Claw is numerical software package for simulating granular-fluid flows, such as landslides, debris flows, and lahars. D-Claw is an extension of the Clawpack codes (www.clawpack.org)
---

---

# Summary

D-Claw is numerical software built from Clawpack codes ([clawpack.org](http://www.clawpack.org)). It is a generalization of GeoClaw ([geoclaw.org](http://www.geoclaw.org)) for solving two-phase flows (granular-fluid mixtures), including landslides, debris flows, and lahars. D-Claw can also be used for simulating hybrid problems that involve the interaction of granular-fluid flows with bodies of water (*eg.*, landslide-generated tsunamis, dam-breach floods, fluid or solid entrainment by inundating flows, such as debris-laden run-off or debris-entraining tsunamis).

The documentation available with Clawpack (v5) and GeoClaw provides a good general overview, however, in addition to unique features, D-Claw is an extension of older versions of Clawpack and GeoClaw (v4.x). (The older v4.x version of Clawpack is available at [github.com/clawpack/clawpack-4.x](https://github.com/clawpack/clawpack4.x). See [clawpack.org](http://www.clawpack.org) for more information.) 

Running D-Claw requires additional set-up parameters. Documentation for D-Claw is currently unsatisfactory, but in progress. Some tips on running it are below. See also the [geoclaw/dclaw-apps](https://github.com/geoclaw/dclaw-apps) on github for examples.


# Obtaining, installing, and running D-Claw

To use D-Claw as is, it is assumed that you have a unix terminal of some kind (*eg.,* linux or Mac OS...for MS Windows you are on your own). 

## Source code

The source code and latest git repository for D-Claw are available on github:

* [github.com/geoflows/D-Claw](https://github.com/geoflows/D-Claw).

A repository for applications is also available:

* [github.com/geoflows/dclaw-apps](https://github.com/geoflows/dclaw-apps).

The application repository is in progress, as is documentation for D-Claw. 

## Installation

"Installation" of D-Claw is essetially just setting some environment variables. Compiling and running D-Claw can then be done using make, with a Makefile in an application directory.  

#### environment variables

The environment variable $CLAW should be set to point to the D-Claw code (the uppermost directory of the source code, "D-Claw/"), or modify the Makefiles if you want to use an environment variable other than $CLAW. For bash shells:
```
export CLAW=/path/to/D-Claw
```

Note: if you are using multiple versions of Clawpack (*eg.,* Clawpack 5.x or GeoClaw and D-Claw), you might want to use, for example, the [environment modules](http://modules.sourceforge.net/) package, which can dynamically set or change your environment under a given shell to ensure that you have a compatible set of paths/versions of software (*eg.*, $PATH, $CLAW, $PYTHONPATH, $MATLABPATH). 

#### python
To use the python set-up scripts that come with D-Claw, you should include $CLAW/python/ in your $PYTHONPATH:
```
export PYTHONPATH=$CLAW/python:$PYTHONPATH
```
You can also then import some tools from the D-Claw python folder if you find them useful:  
```
python> import dclaw
python> import dclaw.topotools as dt 
```
The current version of D-Claw is compatible with python v2.x. If you prefer python version v3.x, see below about developing and contributing to D-Claw.    

#### matlab
If you are going to use matlab plotting with D-Claw, your $MATLABPATH should include $CLAW/matlabgeo, above but in addition to $CLAW/matlab:
```
export MATLABPATH=$CLAW/matlabgeo:$CLAW/matlab:$MATLABPATH
```
For a given application, it is useful to relocate some of the m-files (*eg.,* afterframe.m, setplot2.m, setprob.m, beforeframe.m etc.) included with D-Claw to your local working application directory, where you can modify them to suit you present purposes without modifying your D-Claw source files. 

Note that D-Claw output and these locally modified .m-files (in your application directory) must both be located by matlab, but they are not usually in the same directory. For instance, if the output sub-directory is the current directory where matlab is running, *ie.,*
```
matlab> pwd
/somedir/myapplication/_output
```
then you can issue,
```
matlab> addpath ../
```
to get the local .m-files in the output's parent directory, myapplication, above the same routines in the D-Claw matlab libraries (*ie.* Matlab will add the absolute path for ../ to the top of your path).

## Running D-Claw applications

Running D-Claw applications is very similar to running applications in Clawpack v5.x. See the documentation at [clawpack.org](http://www.clawpack.org). In summary, an application directory houses your Makefile (which correctly resolves the D-Claw source code) and some python set-up scripts. The local file setrun.py in an application folder is used to set runtime parameters, initial conditions and input data (*eg.,* topography DEMs). The make program is used to compile and run the code (as well as produce plots if the python plotting options are used). Plotting can be done with python or matlab scripts, which are included with D-Claw and application examples (*see* *eg.,* [geoclaw/dclaw-apps](https://github.com/geoclaw/dclaw-apps) on github).

Setting up a given simulation essentially amounts to modifying the routine setrun.py. Unfortunately documentation is still lacking for the unique attributes of D-Claw. However, examples that can be modified can be found in the [geoflows/dclaw-apps](https://github.com/geoflows/dclaw-apps) repository.

Further documentation is planned.


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