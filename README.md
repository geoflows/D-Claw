---
title: D-Claw
subtitle: software for granular-fluid flows
description: D-Claw is numerical software package for simulating granular-fluid flows, such as landslides, debris flows, and lahars. D-Claw is an extension of the Clawpack codes (www.clawpack.org)
---

D-Claw
---

D-Claw is numerical software built from Clawpack and Geoclaw codes (see [clawpack.org](http://www.clawpack.org)). D-Claw is a generalization of GeoClaw for solving two-phase flows (granular-fluid mixtures), including landslides, debris flows, and lahars. D-Claw can also be used for simulating hybrid problems that involve the interaction of granular-fluid flows with bodies of water (*eg.*, landslide-generated tsunamis, dam breach floods, fluid or solid entrainment by inundating flows).

The documentation available with Clawpack (v5) and GeoClaw provides a good general overview, however, D-Claw is built on, or extended from, older versions of Clawpack and GeoClaw (v4.x). The older v4.x version of Clawpack is available at [github.com/clawpack/clawpack-4.x](https://github.com/clawpack/clawpack4.x). See [clawpack.org](http://www.clawpack.org) for more information. Running D-Claw requires additional set-up parameters. Documentation for D-Claw is currently unsatisfactory, but in progress.

The source code and latest git repository for D-Claw are available on github:

* [github.com/geoflows/D-Claw](https://github.com/geoflows/D-Claw).

A repository for applications is also available:

* [github.com/geoflows/dclaw-apps](https://github.com/geoflows/dclaw-apps).

The application repository is in progress, as is documentation for D-Claw. If you would like to make contributions to either of these repositories, please follow the development workflow used for Clawpack, described at [www.clawpack.org/developers](http://www.clawpack.org/developers.html). Briefly, please fork the repositories to your own github account, and issue pull requests on a feature branch to github/geoflows, *eg,*:

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

D-Claw inherits the Clawpack licenses and user agreeements. 