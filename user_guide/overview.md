# D-Claw Background and Overview

Katy Barnhart

THIS IS IN PROGRESS. PLACES THAT NEED EXTENSIVE WORK ARE NOTED WITH "TODO".

Lets presume you've correctly installed D-Claw and want to create your own run with your own parameters, input topography, etc. What do you modify and where do you start?

The following summarizes my understanding of the core things that an end-user (e.g., someone who is mostly a simulation-creator, not a source code modifier) would want to know about and change. I came to this understanding by jointly reading:
- Current version of [Clawpack documentation](http://www.clawpack.org/)
- Looking at the D-Claw source code (which reflects my style of learning). I mostly looked at the python interface and plotting scripts.
- Carefully looking through the `dclaw-apps/USGS_Flume/gate_release_example` files and inputs.

Before getting into the details of how one modifies and specifies a D-Claw run, I'll write some background. This is likely most usesful to someone who has never used anything in the Clawpack ecosystem before.

## 1. General Clawpack, Geoclaw, D-Claw introduction
The basic D-Claw workflow inherits from Clawpack workflows. There are four main clawpack solvers:

- Classic, the initial fortran from RJlV book
- AMRClaw, adaptive mesh refinement in fortran
- GeoClaw, builds on AMRClaw to add some special algorithms for geophysics.
- PyClaw includes some additional new algorithms.

Parallelism is possible in all. Fortran uses OpenMP and PyClaw uses MPI and PETSc.

The [current documentation of Clawpack](http://www.clawpack.org/contents.html) is very comprehensive.

[This page]() goes to the Clawpack v4.X.X, on which D-Claw currently rests.

### 1.1 Very brief summary of D-Claw

The two key papers that describe the derivation of the D-Claw equations and their numerical solution:

Iverson, R. M., & George, D. L. (2014). A depth-averaged debris-flow model that includes the effects of evolving dilatancy. I. Physical basis. Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences, 470(2170), 20130819. https://doi.org/10.1098/rspa.2013.0819

George, D. L., & Iverson, R. M. (2014). A depth-averaged debris-flow model that includes the effects of evolving dilatancy. II. Numerical predictions and experimental tests. Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences, 470(2170), 20130820. https://doi.org/10.1098/rspa.2013.0820

### 1.2 What are the state variables?
Five dependent variables (see George and Iverson (2014) Section 2a):
- h(x, y, t)
- hu(x, y, t)
- hv(x, y, t)
- hm(x, y, t)
- pb(x, y, t)

where h is the height normal to the bed of a virtual free surface, u and v are depth averaged flow velocities in the x and y directions, m is the depth-averaged solid volume fraction, and pb is the pore fluid pressure at the basal boundary.

The bed (b) is additionally specified, and (optionally) an erodible layer of thickness (a)
### 1.3 Conceptual description of inputs:

In order to run, D-Claw needs the following information:
- Specification of the initial topographic surface (eta), thickness of mobile material (h), and bed (b). These three quantities are related (eta = b + h) so only two need to be specified.
- Optionally a thickness (a) of erodible material can be specified.
- Optionally a slope-normal (theta) can be specified.
- A number of values scalar values (described in the table below) related to physical properties. Some of these may be specified as scalars or spatially variable values.

In addition, a number of standard Clawpack inputs are required. The most core being:
- Definition of the computational grid, including AMR information.
- Definition of the simulation timeframe.
- Statement of when, where, and at what resolution the simulation state should be written to an output file.  
- Definition of any gauges (see [this Clawpack page](http://www.clawpack.org/gauges.html)).

Note that some inputs must be specified as spatially distributed (e.g., h, eta) while some others (e.g., initial value for m, coulomb normal friction angle (phi)) may be specified as either a scalar or spatially distributed.

Note also that the way D-Claw controls AMR refinement is different than how Clawpack handles it. See the discussion of "flowgrades" below.

### 1.4 Conceptual statement of failure options

TODO.

### 1.5 The `q` and `aux` arrays

There are two core arrays that Clawpack maintains. One is called `aux` and refers to non time-variable quantities. One is called `q` and it refers to the state variables. `q` often also contains other time-variable quantities that are commonly desired for writing output.

When Clawpack refers to writing the outputs of q (e.g., at a gauge) this refers to:

q = (h, hu, hv, hm, pb)

Note, however, that in D-Claw the output written in the fort.qXXX files includes more variables than these five state variables. There are eight variables written to the q vector and they are as follows:
- h
- hu
- hv
- hm
- pb
- something related to the optional particle size distribution function which permits treatment of segregation.
- positive downward elevation change of the erodible layer, a
- eta, the surface elevation.

The aux array contains the following:
- topography, b
- elements two and three of the aux array relate to how geoclaw handles latitude and longitude.
- see element two.
- phi, the coulomb friction angle.
- theta, the bed normal slope angle (if specified) OR a, the thickness of erodible material. Only one of these two options can be used at a time.
- fs (this and the following three elements of aux are related to the factor of safety used in the Riemann solver).
- fsphi
- taudir_x
- taudir_y

### 1.6 A note about coordinate system and slope normals

George and Iverson (2014) Section 2c provide a discussion about the numerical challenges that come from using a depth-averaged model with complex basal topography. The USGS flume example shows how to specify theta using the aux array to use a slope normal coordinate system. This, however is only possible with such simple topography. This section exists to remind users that this is an important issue and that this section of the paper provides some background.

### 1.7 Grain size segregation.
Optional, TODO.

### 1.8 Best practices with respect to topography, grid definition, and AMR

TODO: Want to summarize here best practices for how desired/anticipated AMR, grid definition, and topography input scales all interact.

## 2. Specification of Inputs and Running.

TODO: Finish.

This is done in `setrun.py` and `setplot.py`. Run using make .output. See info in Clawpack docs about Makefile approach.

`setrun.py` itself has lots of notes in it (some of which are copied into here) I've tried to synthesize the core bits and discuss the ones that were not immediately self explanatory (for me). The order and grouping here is also distinct from that in `setrun.py`.

### 2.1 `setrun.py` used to set parameters

- Clawpack docs on specifying these parameters in [classic](http://www.clawpack.org/setrun.html#setrun) and [geoclaw](http://www.clawpack.org/setrun_geoclaw.html#setrun-geoclaw). **Start by reading these**
- There is a nice distinction between what parameters are specific to clawpack (`clawdata.attribute=value`), geoclaw (`geodata.attribute=value`), and D-claw (`digdata.attribute=value`)
- In addition to the clawpack and geoclaw parameters, DCLAW can take the following additional parameters. The full list of these and their defaults can be found in the definition of the `DigclawInputData` (line 1427 of  `D-Claw/python/pyclaw/data.py`).

### 2.2 D-Claw specific scalars

Many of these are defined in George and Iverson (2014). Where relevant there is a place for this definition to be noted (TODO).
Note that some of these can be set as spatially variable using values of `q` or `aux`. This is accomplished setting auxinit or qinit (see below). Here q_1 means the first element of the q array (python index zero, iqinit=1)

| Parameter Name     | How to specify if spatially variable? | Default Value | Description (units)| Where to find in G&I(2014) |
| --- | --- | --- | --- | --- |
| rho_s            |         | 2700.0        | solid grain density (kg/m^3) |    |
| rho_f            |         | 1000.0        | pore-fluid density  (kg/m^3) |    |
| phi_bed          | aux_4   | 40.           | basal friction angle (degrees) |    |
| theta_input      | aux_5   | 0.            | slope angle (degrees) |    |
| delta            |         | 0.01          | characteristic grain diameter (m) |    |
| kappita          |         | 0.0001        | characteristic grain diameter parameter in kperm (m) |    |
| mu               |         | 0.001         | viscosity of pore-fluid (Pa-s) |    |
| alpha_c          |         | 1.0           | debris compressibility constant (#) |    |
| m_crit           |         | 0.62          | critical state value of m (#) |    |
| c1               |         | 1.0           | dilation coefficient 1 (#) |    |
| m0               | q_5     | 0.52          | initial solid volume fraction (#) |    |
| sigma_0          |         | 1.e3          | baseline stress for definition of compessibility |    |
| alpha_seg        |         | 0.0           | coefficient of segregation velocity profile |    |
| init_ptype       |         | 0             | 0 = hydrostatic; 1 or 2 = failure pressure; 3 or 4 = rising pressure |    |
| init_pmax_ratio  |         | 1.0           | pressure rise to hydrostatic times init_pmax_ratio |    |
| init_ptf         |         | 1.0           | pressure will rise until t = init_ptf without dilatancy |    |
| init_ptf2        |         | 0.0           | will rise until t = init_ptf2 |    |
| bed_normal       |         | 0             | bed_normal = 1 requires theta in aux. for slope in one direction |    |
| phi_seg_coeff    |         | 0.0           | adjustment to friction coefficient based on segregation |    |
| entrainment      |         | 0             | flag for entrainment, 0 = no entrainment, specify thickness of a using aux_5 |    |
| entrainment_rate |         | 0.2           | rate of entrainment parameter 0-1 |    |

These are specified by modifying the object `digdata` in `setrun.py`, for example, to set the pore fluid density to 1100 one might do:
```python
digdata.rho_f = 1100.0
```

### 2.3 Format of spatially variable inputs
Start by reading the information in the [topography data documentation on Clawpack website](http://www.clawpack.org/topo.html#topo). Cribnotes: Four styles of input are supported (3 ascii and 1 netcdf).
- There are multiple types of topography-like files but at core they are all (x,y,z). Esri ascii format is topotype 3.
- Another good page is the [geoclaw doc page on topography data file parameters](http://www.clawpack.org/setrun_geoclaw.html#topography-data-file-parameters) for definition of each input.
- Clawpack has a lot of python tools for working with topography-like data. Start [here](http://www.clawpack.org/topotools_module.html)
In this document, ee also the section on grid definition and how grid resolution and topographic resolution are related (Section 2.5).


TODO: Make clear that there are some indexing differences between the current clawpack docs and DCLAW. So there may be places where the Clawpack docs are not quite right.

### 2.4 Different types of spatially variable inputs that might be set.

#### 2.4.1 Topography (b)

Topography specified in `geodata.topofiles`. [Clawpack doc link](http://www.clawpack.org/setrun_geoclaw.html#topography-data-file-parameters) which includes definitions of parameters:
```python
geodata.topofiles = []
geodata.topofiles.append([topotype, minlevel, maxlevel, t1, t2, fname])
```
in short
- `topotype`, topography format type
- `minlevel`, minimum refinement level enforced in this region between times `t1` and `t2`.
- `maxlevel `, maximum refinement level enforced in this region between times `t1` and `t2`.
- `t1`, start time of refinement control time period.
- `t2`, end time of the refinement control time period.
- `fname`, filename.

In addition, topography can be specified at multiple scales. From the docs:

> More than one topo file can be specified (see Topography data file parameters) that might cover overlapping regions at different resolutions. The union of all the topo files should cover the full computational domain specified (and may extend outside it). Internally in GeoClaw a single piecewise-bilinear function is constructed from the union of the topo files, using the best information available in regions of overlap. This function is then integrated over computational grid cells to obtain the single topo value in each grid cell needed when solving depth averaged equations such as the shallow water equations with these finite volume methods. Note that this has the feature that if a grid cell is refined at some stage in the computation, the topo used in the fine cells have an average value that is equal to the coarse cell value. This is crucial in maintaining the ocean-at-rest steady state, for example.

#### 2.4.2. Displacement
This is unlikely to be used by D-Claw and is meant more for tsunami initiation.

```python
geodata.dtopofiles = []
geodata.topofiles.append([dtopotype, minlevel, maxlevel, fname])
```
See [this page](http://www.clawpack.org/topo.html#topo-dtopo) for detailed info on the dtopotypes.

#### 2.4.3 qinit
The qinit specification is used to state modification of the initial conditions represented in the q array. [Geoclaw doc description link](http://www.clawpack.org/setrun_geoclaw.html#setrun-qinit). These are specified as:

```python
geodata.qinitfiles = []
geodata.qinitfiles.append(qinitftype,iqinit, minlevel, maxlevel, fname])
```
The `qinitftype` is the file type, the same as is permitted for topo. `iqinit` is the type of perterbation, e.g. which element of q is perterbed. Note here that while h is the first element of q, with python index zero, it cooresponds to `iqinit=1`

While the surface elevation, eta, is not technically a state variable, it can be set as the eight element of qinit.

TODO: what happens if ove- specified (e.g., topo, h, and eta), or `digdata.m0` and a value for q_5.
TODO: related. If all others specified, then `p_b` implied, yes?

#### 2.4.4 auxinit
This is where you set initial values for the aux array. Same as qinit.

```python
geodata.auxinitfiles = []
geodata.auxinitfiles.append([auxinitftype, iauxinit, minlevel, maxlevel, fname])
```

### 2.5 Definition of the grid scale and AMR levels.

TODO.

### 2.6 Definition of simulation time.

TODO

### 2.7 Definition of output writing.
TODO. Include Turn on/off aux?


### 2.8 Flowgrades and AMR control
TODO.
This controls AMR refinement.

`geodata.flowgrades` from the input file:
  ```python
  # for using flowgrades for refinement append lines of the form
  # [flowgradevalue, flowgradevariable, flowgradetype, flowgrademinlevel]
  # where:
  #flowgradevalue: floating point relevant flowgrade value for following measure:
  #flowgradevariable: 1=depth, 2= momentum, 3 = sign(depth)*(depth+topo) (0 at sealevel or dry land).
  #flowgradetype: 1 = norm(flowgradevariable), 2 = norm(grad(flowgradevariable))
  #flowgrademinlevel: refine to at least this level if flowgradevalue is exceeded.
  geodata.flowgrades = []
  geodata.flowgrades.append([flowgradevalue, flowgradevariable, flowgradetype, flowgrademinlevel])
  ```
### 2.9 Source terms and  the source  .f90 file

TODO

## 3. Output and postprocessing with `setplot.py`
Mostly various textfiles.

TODO

### 3.1 Output files and format.
TODO
### 3.2 `setplot.py`
TODO. Prob just point to Clawpack docs as this really is a language.

## 4. Organization of the repository

In the D-Claw repository we have the following key directories:
- Fortran source
  - `amrclaw` (adaptive mesh refinement claw)
  - `clawpack` (core claw, this is v4.X.X of Clawpack)
  - `geoclaw`
      - The core of D-Claw is located in `D-Claw/geoclaw/2d/lib_dig`

- `python`, which contains two unpackaged python modules and some additional tools. These are made accessible to the runs by setting the `PYTHONPATH` to include them.

  - `dclaw`. Upon initial review, this is mostly code to deal with D-Claw specific pre or postprocessing. It has quite a few comments that are helpful to understand the D-Claw inputs and outputs.

  - `pyclaw` The python wrapping of clawpack. This provides the ability to specify and run clawpack solvers. The modern (clawpack v5.6) version has a packaged version of pyclaw, compatible with python 3.6.
