---
title: D-Claw
subtitle: software for granular-fluid flows
description: D-Claw is numerical software package for simulating granular-fluid flows, such as landslides, debris flows, and lahars. D-Claw is an extension of the Clawpack codes (www.clawpack.org)
---

D-Claw
---

# Summary

D-Claw is numerical software built from Clawpack and Geoclaw codes (see [clawpack.org](http://www.clawpack.org)). D-Claw is a generalization of GeoClaw for solving two-phase flows (granular-fluid mixtures), including landslides, debris flows, and lahars. D-Claw can also be used for simulating hybrid problems that involve the interaction of granular-fluid flows with bodies of water (*eg.*, landslide-generated tsunamis, dam breach floods, fluid or solid entrainment by inundating flows).

The documentation available with Clawpack (v5) and GeoClaw provides a good general overview, however, D-Claw is built on, or extended from, older versions of Clawpack and GeoClaw (v4.x). The older v4.x version of Clawpack is available at [github.com/clawpack/clawpack-4.x](https://github.com/clawpack/clawpack4.x). See [clawpack.org](http://www.clawpack.org) for more information. Running D-Claw requires additional set-up parameters. Documentation for D-Claw is currently unsatisfactory, but in progress.


# Use

## source code

The source code and latest git repository for D-Claw are available on github:

* [github.com/geoflows/D-Claw](https://github.com/geoflows/D-Claw).

A repository for applications is also available:

* [github.com/geoflows/dclaw-apps](https://github.com/geoflows/dclaw-apps).

The application repository is in progress, as is documentation for D-Claw. 

## Tips for running D-Claw

##### environment variables

* make sure the environment variable $CLAW points to the D-Claw code (the uppermost directory of the repository, "D-Claw.") (Or modify the Makefiles if you want to use an environment variable other than $CLAW. For bash shells:
```
export CLAW=/somedir/D-Claw
```

* if you are using multiple versions of Clawpack (*eg.,* Clawpack 5.x or GeoClaw and D-Claw), you might want to use, if you don't already, the [environment modules](http://modules.sourceforge.net/) package, which can dynamically set or change your environment under a given shell, to make sure you have a compatible set of paths/versions of software (*eg.*, $PATH, $CLAW, $PYTHONPATH, $MATLABPATH).

##### python
* its a good idea to have $CLAW/python/ in your $PYTHONPATH:
```
export PYTHONPATH=$CLAW/python:$PYTHONPATH
```
You can then import some tools from the D-Claw python folder if they are useful:  
```
python> import dclaw
python> import dclaw.topotools as dt 
```  

##### matlab
* if you are going to use matlab plotting with D-Claw, be sure that your $MATLABPATH contains $CLAW/matlabgeo and then $CLAW/matlab
```
export MATLABPATH=$CLAW/matlabgeo:$CLAW/matlab:$MATLABPATH
```
* if you are using matlab and "plotclaw2.m", note that output and local .m-files must both be located by matlab, but they are not usually in the same directory. If the output directory is your current directory where matlab is running, you can
```
matlab> addpath ../
```
to get the local .m-files in the output's parent directory correctly on your path. (Matlab will add the absolute path for ../)


* if you are using multiple versions of Clawpack (*eg.,* Clawpack 5.x or GeoClaw and D-Claw), you might want to use, if you don't already, the [environment modules](http://modules.sourceforge.net/) package, which can dynamically change your environment variables (such as $PATH, $CLAW, $PYTHONPATH, $MATLABPATH) for your running shell, to make sure all paths are correct and compatible.


## development

If you would like to make contributions to D-Claw or dclaw-apps, please follow the development workflow used for Clawpack, described at [www.clawpack.org/developers](http://www.clawpack.org/developers.html). Briefly, please fork the repositories to your own github account, and issue pull requests on a feature branch to github/geoflows, *eg,*:

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