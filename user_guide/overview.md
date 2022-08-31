# D-Claw Background and Overview

Katy Barnhart

THIS IS IN PROGRESS. PLACES THAT NEED EXTENSIVE WORK ARE NOTED WITH "TODO".

Lets presume you've correctly installed D-Claw and want to create your own run with your own parameters, input topography, etc. What do you modify and where do you start?

The following summarizes my understanding of the core things that an end-user (e.g., someone who is mostly a simulation-creator, not a source code modifier) would want to know about and change. I came to this understanding by jointly reading:
- Current version of [Clawpack documentation](http://www.clawpack.org/)
- Looking at the D-Claw source code (which reflects my style of learning). I mostly looked at the python interface and plotting scripts.
- Carefully looking through the `dclaw-apps/USGS_Flume/gate_release_example` files and inputs.
- Asking Dave George lots of questions, which he thoughtfully answered.

Before getting into the details of how one modifies and specifies a D-Claw run, I'll write some background. This is likely most useful to someone who has never used anything in the Clawpack ecosystem before (which is where I started).

## 1. General Clawpack, Geoclaw, D-Claw introduction
The basic D-Claw workflow inherits from Clawpack workflows. There are four main clawpack solvers:

- Classic, the initial fortran from RJlV book
- AMRClaw, adaptive mesh refinement in fortran
- GeoClaw, builds on AMRClaw to add some special algorithms for geophysics.
- PyClaw includes some additional new algorithms.

Parallelism is possible in all. Fortran uses OpenMP and PyClaw uses MPI and PETSc. I usually operate embarrassingly parallel efforts because I like to explore many parameter combinations. Thus I've not engaged with any parallelization.

The [current documentation of Clawpack](http://www.clawpack.org/contents.html) is very comprehensive.

[This page](https://depts.washington.edu/clawpack/users-4.6/index.html) goes to the Clawpack v4.6, on which D-Claw approximately rests. TODO, correct version?

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
- Optionally a slope-normal (theta) can be specified. Note that only a or theta can be specified at the same time given array allocation set up (they occupy the same place in the aux array).
- A number of values scalar values (described in the table below) related to physical properties. Some of these may be specified as scalars or spatially variable values.

In addition, a number of standard Clawpack inputs are required. The most core being:
- Definition of the computational grid, including AMR information.
- Definition of the simulation timeframe.
- Statement of when, where, and at what resolution the simulation state should be written to an output file.  
- Definition of any gauges (see [this Clawpack page](http://www.clawpack.org/gauges.html)).

Note that some inputs must be specified as spatially distributed (e.g., h, eta) while some others (e.g., initial value for m, coulomb normal friction angle (phi)) may be specified as either a scalar or spatially distributed.

Note also that the way D-Claw controls AMR refinement is different than how Clawpack handles it. See the discussion of "flowgrades" below.

### 1.4 The `q` and `aux` arrays

There are two core arrays that Clawpack maintains. One is called `aux` and refers to non time-variable quantities. One is called `q` and it refers to the state variables. `q` often also contains other time-variable quantities that are commonly desired for writing output.

When Clawpack refers to writing the outputs of q (e.g., at a gauge) this refers to:

q = (h, hu, hv, hm, pb)

Note, however, that in D-Claw the output written in the fort.qXXX files includes more variables than these five state variables. There are eight variables written to the q vector and they are as follows:
- h
- hu
- hv
- hm
- pb
- h*chi, where chi is the fraction of particle species 1. This is optionally used in the treatment of particle species segregation.
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

### 1.5 Initialization of basal pressure pb (q5)

Basal pressure, pb (q5) is initialized in the subroutine qinit(...) located in the file geoclaw/2d/lib_dig/qinit_geo.f.

The method for setting basal pressure is controlled by the variable `init_ptype` variable. Depending on the value for `init_ptype` some other variables (`init_p***` variables may or may not be used).

#### 1.5.0 How is init_pmin_ratio calculated

The parameter `init_pmin_ratio` is fraction of hydrostatic pressure used to initialize spatially variable pressure if `init_ptype>0`. It depends on the gradient of the surface (eta) and the basal friction angle (phi). It does not depend on m because later on (todo: in what routine) asking whether there is sufficent m to keep material static is addressed.

Across the entire domain the ratio between surface slope and tan(phi) is calculated at all locations where h>drytolerance.

If `init_ptype` is 1 or 3, `init_pmin_ratio` is the absolute minimum pressure ratio that yields liquifaction at at least one location in the domain. That is, if multiple cells may fail at the same lowest pressure ratio, then all of those cells, but no other cells will fail.

If type 2 or 4, `init_pmin_ratio` is the average value of grad_eta/tan(phi) across the domain.

grad_eta = sqrt(d_eta_dx^2 + d_eta_dy^2)
init_ptype_ratio = 1.0 - (grad_eta / tan (phi))

Main options:

#### 1.5.1 Zero pressure or user defined files.
If `init_ptype = -1` then q5 is initialized to zero everwhere. If files for pb are also provided in the qinitfiles list, then the value of pb=0 is overwritten by the quantity specified times q1. That is, the value provided is the pressure relative to material thickness pb/h.

Note, if user provided files are provided but `init_ptype != -1` then the user provided files will be overwritten without any warning.

#### 1.5.2 Hydrostatic
If `init_ptype = 0` then the basal pressure is set to hydrostatic everwhere based on the material thickness h and the fluid pressure rho_f:

pb = rho_f * gmod * h

Note that in this calculation gmod is gravitational acceleration normal to the basal surface.

#### 1.5.3 Set to Failure Pressure
If `init_ptype = 1 or 2` then the

if h > drytolerance:
  rho = sv*rho_s + (1.0-sv)*rho_f
  pb = init_pmin_ratio * rho * gmod * h

1 = minimum <- gives you failure at t=0 at the one location in which the lowest pressure yields liquifaction.

2 = average <- gives you failure at t=0 at a range of locations in which the pressure needed for liquifaction is lower than average failure pressure.

#### 1.5.4 Rising Pressure
`init_ptype = 3 or 4`

The pressure is raised linearly in time over the duration t=0 to t=init_ptf from pressure = zero to pressure = pfail where

pfail = hydrostatic * init_pmin_ratio

This is acomplished in the routine b4step2(...) (geoclaw/2d/lib_dig/b4step2_geo.f).

Like in the failure pressure options the value for init_pmin_ratio can be set based on either the minimum or average value.

3 = minimum
4 = average

Setting init_ptype = 3 should be equivalent init_ptype = 1, but with the time of failure at t=init_ptf instead of t=0 (this has not been verified).

In contrast, type 2 and 4 should be different because in init_ptype = 2 many locations may fail at t=0, while in init_ptype=4, failure may occur at different locations over the duration t=0 to t=init_ptf. Note (2022-08-09), as best as KRB/RPJ can tell from the code itself, use of type=4 would result in resetting the pressure back to the linear increase *even in cells which have failed and are moving* this could be fixed by adding a velocity check in b4step2. The use of init_ptype = 4 is not recommended.

#### 1.5.4 Misc notes:
NOTES from DG: For future development the (non-qinitfile) options should be 0, hydrostatic, failure (which is P_scaled = min_domain (p_b/h) such that failure occurs , applied everywhere: p_b = P_scaled*h). All of these options could occur at t=0 with no time interval used...Qinit files should override every other option.  

2022-08-09 : KRB cannot find evidence that any option uses `init_pmax_ratio`. It is read, but is only reference by geoclaw/1d/lib_dig/b4step1_dig.f. This routine not used in current makefiles.

### 1.6 A note about coordinate system and slope normals

George and Iverson (2014) Section 2c provide a discussion about the numerical challenges that come from using a depth-averaged model with complex basal topography. The USGS flume example shows how to specify theta using the aux array to use a slope normal coordinate system. This, however is only possible with such simple topography. This section exists to remind users that this is an important issue and that this section of the paper provides some background.

### 1.7 Grain size segregation.
Optional, TODO.

### 1.8 Best practices with respect to topography, grid definition, and AMR

TODO: Want to summarize here best practices for how desired/anticipated AMR, grid definition, and topography input scales all interact.

Key concept is that the AMR gridding is topography resolution agnostic. If you define topography at a resolution of X and amr refines to 0.1X it will interpolate based on the underlying spline methods. If you define topo at X and AMR is operating at 10X it will conservatively represent an integral of the topo at its resolution.

Only way in which AMR "cares" about topo resolution is that if multiple topography files are provided it will use the highest resolution one at any given space.

Note, the min and max AMR levels defined when topo and qinit files are initialized define AMR "regions" of the permissible refinement levels. See discussion of "regions" below.

If more than one topography file is provided at the same resolution for an overlapping part of the domain, the file that is provided last to the list will be used.

### 1.8 Overspecification

Text from Dave: Redundancy/overspecification: there is no error catching for this, which might be useful in the future. However, it would require location overlap checks, because often a variable might be set by different methods for different parts of the domain. (e.g., it might be useful to set water depth for a lake using eta and landslide material using h, etc.). Ideally, qinitfiles should override the default values...however, pressure will need some attention to make sure that's the case...currently, it should be the case for h,hu,hv,m. There are no error checks for incompatible qinitfiles...that could be useful and not too painful to do with qinitfile overlap checks. Currently, in regions where there are no qinitfiles, the defaults are used (aside from pressure, which is about to undergo changes to handle the very issues you raise).

## 2. Specification of Inputs and Running.

TODO: Finish.

This is done in `setrun.py` and `setplot.py`. Run using make .output. See info in Clawpack docs about Makefile approach.

`setrun.py` itself has lots of notes in it (some of which are copied into here) I've tried to synthesize the core bits and discuss the ones that were not immediately self explanatory (for me). The order and grouping here is also distinct from that in `setrun.py`.

### 2.1 `setrun.py` used to set parameters

- Clawpack docs on specifying these parameters in [classic](http://www.clawpack.org/setrun.html#setrun) and [geoclaw](http://www.clawpack.org/setrun_geoclaw.html#setrun-geoclaw). **Start by reading these**
- There is a nice distinction between what parameters are specific to clawpack (`clawdata.attribute=value`), geoclaw (`geodata.attribute=value`), and D-claw (`digdata.attribute=value`)
- In addition to the clawpack and geoclaw parameters, DCLAW can take the following additional parameters. The full list of these and their defaults can be found in the definition of the `DigclawInputData` (line 1427 of  `D-Claw/python/pyclaw/data.py`).

### 2.2 Default values for spatially variable inputs (q and aux)

depth is sea level (unless specified by qinitfiles)
both momentum terms set as zero.

pressure not typically set as spatially variable.
m0 is typically a scalar.

### 2.3 D-Claw specific scalars

Many of these are defined in Iverson and George (2014) and George and Iverson (2014). Where relevant there is a place for this definition to be noted (TODO).
Note that some of these can be set as spatially variable using values of `q` or `aux`. This is accomplished setting auxinit or qinit (see below). Here q_1 means the first element of the q array (python index zero, iqinit=1)


For the most up-to-date description of D-Claw specific inputs, look at the dictionary _DIG_ATTRS python/pyclaw/data.py


| Parameter Name     | How to specify if spatially variable? | Default Value | Description (units)| Where to find in G&I(2014) |
| --- | --- | --- | --- | --- |
| rho_s            |         | 2700.0        | solid grain density (kg/m^3) |    |
| rho_f            |         | 1000.0        | pore-fluid density  (kg/m^3) |    |
| phi_bed          | aux_4   | 40.           | basal friction angle (degrees) |    |
| theta_input      | aux_5   | 0.            | slope angle (degrees) |    |
| delta            |         | 0.01          | characteristic grain diameter (m) |    |
| kappita          |         | 0.0001        | permeability at m=setdig.m0 (m^2) | k0 in G&I eq 2.7 - modified such that setdig.m0 is used instead of 0.6  |
| mu               |         | 0.001         | viscosity of pore-fluid (Pa-s) |    |
| alpha_c          |         | 1.0           | debris compressibility constant (#) | This is the a constant in the equation for alpha  |
| m_crit           |         | 0.62          | critical state value of m (#) |    |
| c1               |         | 1.0           | dilation coefficient 1 (#) | This is a tuning nob. Not in an equation. Leave at 1   |
| m0               | q_5     | 0.52          | initial solid volume fraction (#) |    |
| sigma_0          |         | 1.e3          | baseline stress for definition of compressibility | Defined in G/I eq 2.8 |
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

### 2.4 Format of spatially variable inputs
Start by reading the information in the [topography data documentation on Clawpack website](http://www.clawpack.org/topo.html#topo). Cribnotes: Four styles of input are supported (3 ascii and 1 netcdf).
- There are multiple types of topography-like files but at core they are all (x,y,z). Esri ascii format is topotype 3.
- Another good page is the [geoclaw doc page on topography data file parameters](http://www.clawpack.org/setrun_geoclaw.html#topography-data-file-parameters) for definition of each input.
- Clawpack has a lot of python tools for working with topography-like data. Start [here](http://www.clawpack.org/topotools_module.html)
In this document, ee also the section on grid definition and how grid resolution and topographic resolution are related (Section 2.5).

TODO: Make clear that there are some indexing differences between the current clawpack docs and DCLAW. So there may be places where the Clawpack docs are not quite right.

### 2.5 Different types of spatially variable inputs that might be set.

#### 2.5.1 Topography (b)

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

Note. There may be some line ending requirements for .tt3 files. I'm not 100% sure, but I've had some issues fortran reading in topo files with files made on windows that are "fixed" when I read a .tt3 file into a python array and then back out to .tt3 using the topotools.

#### 2.5.2. Displacement
This is unlikely to be used by D-Claw and is meant more for tsunami initiation.

```python
geodata.dtopofiles = []
geodata.topofiles.append([dtopotype, minlevel, maxlevel, fname])
```
See [this page](http://www.clawpack.org/topo.html#topo-dtopo) for detailed info on the dtopotypes.

#### 2.5.3 qinit
The qinit specification is used to state modification of the initial conditions represented in the q array. [Geoclaw doc description link](http://www.clawpack.org/setrun_geoclaw.html#setrun-qinit). These are specified as:

```python
geodata.qinitfiles = []
geodata.qinitfiles.append(qinitftype,iqinit, minlevel, maxlevel, fname])
```
The `qinitftype` is the file type, the same as is permitted for topo. `iqinit` is the type of perturbation, e.g. which element of q is perterbed. Note here that while h is the first element of q, with python index zero, it cooresponds to `iqinit=1`

While the surface elevation, eta, is not technically a state variable, it can be set as the eighth element of qinit.

NOTE: q4 must be provided as the quantity m, not hm, and q5 must be provided as the quantity pb/h, not pb.

TODO: What happens if inputs are over specified in a conflicting way (e.g., topo, h, and eta such that topo + h != eta), or `digdata.m0` and a value for q_5.

#### 2.5.4 auxinit
This is where you set initial values for the aux array. Same as qinit.

```python
geodata.auxinitfiles = []
geodata.auxinitfiles.append([auxinitftype, iauxinit, minlevel, maxlevel, fname])
```

### 2.6 Definition of the grid scale and AMR Levels

Domain extent and number of grid cells are set as follows:
```python
# Lower and upper edge of computational domain:
clawdata.xlower = -10.0
clawdata.xupper =  140.0
clawdata.ylower =  -6.0
clawdata.yupper =   8.0


# Number of grid cells:
clawdata.mx = 200
clawdata.my = 28
```

AMR levels given as:
```python
# ---------------
# AMR parameters:
# ---------------

# max number of refinement levels:
mxnest = 2

clawdata.mxnest = -mxnest   # negative ==> anisotropic refinement in x,y,t

# List of refinement ratios at each level (length at least mxnest-1)
clawdata.inratx = [5,5,2,4]
clawdata.inraty = [5,5,2,4]
clawdata.inratt = [5,5,2,4]
```

AMR levels are 1 indexed. AMR level 1 is coarsest.

If variable t refinement is used, then inratt is only used to help improve the "guess" at time refinement when moving from coarse to fine grids. If you don't have a really good reason to make these three inratX values different, you should keep them the same.

[Link to example AMRCLAW setrun.py from v4.6](http://depts.washington.edu/clawpack/users-4.6/amrclaw/setrun_amrclaw_sample.html?highlight=mxnest)

### 2.7 Definition of simulation time.

Initial time is set with:
```python
clawdata.t0 = 0.0
```

A few ways to set final time, it is implied by the output times and will stop once the final output time has been reached.

```python
clawdata.tfinal = 20.0
```

Time stepping and courant conditions  described by:
```python
# --------------
# Time stepping:
# --------------

# if dt_variable==1: variable time steps used based on cfl_desired,
# if dt_variable==0: fixed time steps dt = dt_initial will always be used.
clawdata.dt_variable = 1

# Initial time step for variable dt.
# If dt_variable==0 then dt=dt_initial for all steps:
clawdata.dt_initial = 1.e-16

# Max time step to be allowed if variable dt used:
clawdata.dt_max = 1e+99

# Desired Courant number if variable dt used, and max to allow without
# retaking step with a smaller dt:
clawdata.cfl_desired = 0.25
clawdata.cfl_max = 0.9

# Maximum number of time steps to allow between output times:
clawdata.max_steps = 100000
```

Note, in order for variable dt to be used the geoclaw flag `geodata.variable_dt_refinement_ratios` must be set to `True`. This `clawdata.dt_variable` is an older flag that is trumped by the variable dt refinement.

Also, if variable dt refinement is used the dt_initial seems to be ignored. That is an initial value for dt at the highest level is identified using the CFL condition and the wavespeeds rather than this dt_initial.

Also, the dt_max refers to the max at the coarsest AMR level. It not clear if its presently being honored, but we are working on it over at [PR 13](https://github.com/geoflows/D-Claw/pull/13).

The reason for two Courant numbers is written right above the variables in setrun.py. There is computational expense in both taking small timesteps and in taking a timestep that is too big such that the wavespeeds known at the end of the timestep mean that the timestep was in violation of the cfl_max value. So giving a range between cfl_desired and cfl_max means that D-Claw will aim for clf_desired and won't need to retake timesteps so long as it doesn't exceed cfl_max. Remember, end of timestep wavespeeds are not know at the beinning of the timestep.

### 2.8 Definition of output writing.

In the Clawpack v5.7 docs the following options are listed in the [example AMRClaw setrun.py file](http://www.clawpack.org/setrun_amrclaw_sample.html#setrun-amrclaw-sample).

```python
clawdata.output_format = 'ascii'       # 'ascii', 'binary', 'netcdf'

clawdata.output_q_components = 'all'   # could be list such as [True,True]
clawdata.output_aux_components = 'none'  # could be list
clawdata.output_aux_onlyonce = True    # output aux arrays only at t0
```

### 2.9 Flowgrades and AMR control

There are two ways that AMR refinement is controlled in D-Claw. The first is a set of AMR regions which define the permissible range of AMR refinement over a region. Note that if you have a large region which permits refinement between levels 1 and 4, and then a smaller contained region which indicates refinement between levels 2 and 3, the larger,  more permissive region with trump. In other words, regions can't be used to "restrict" AMR levels, unless they are the ONLY region that applies to an area.

Regions are specified as bounding boxes.

```python
geodata.regions=[]
geodata.regions.append([minlevel,maxlevel,t1,t2,x1,x2,y1,y2])
```
where
`minlevel` and `maxlevel` are the minimum and maximum AMR levels.
`t1` and `t2` are the bounding times
and `x1`, `x2`, `y1`, `y2` define the bounding box.

Providing a topography or qinit file implicitly defines a region around these data sources.

In addition to regions D-Claw controls refinement through flowgrades. Refinement occurs where one of the three flow grade variables exceeds a specified flow grade value.

Specified as
```python
geodata.flowgrades = []
geodata.flowgrades.append([flowgradevalue, flowgradevariable, flowgradetype, flowgrademinlevel])
```
where
`flowgradevalue` is the floating point relevant flowgrade value for following measure
`flowgradevariable` is one of the following
    - 1 = depth
    - 2 = momentum,
    - 3 = `sign(depth)*(depth+topo)` (0 at sealevel or dry land).
`flowgradetype` is one of the following:
    - 1 = norm(flowgradevariable)
    - 2 = norm(grad(flowgradevariable))

If the `flowgradevalue` is exceeded then D-Claw will refine up to at least the level indicated by `flowgrademinlevel` unless refinement is limited by AMR refinement limits defined by AMR regions.

If no flowgrades are specified, then D-Claw will use "tsunami-style" refinement based on the values specified in `geodata.depthdeep` and `geodata.maxleveldeep`. These define an area adjacent to the shoreline (up to waterlevel of `geodata.depthdeep`) where refinement can occur as specified by the regions. In deeper water, the level is limited to `geodata.maxleveldeep`.

Flowgrades and tsunami-style refinement are implemented in a disjoint way. That is, you can use either flowgrades OR tsunami-style refinement, but not both. If flowgrades are specified, they will take precedence and tsunami-style refinement parameters will be ignored.

### 2.10 Source terms and  the source .f90 file

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
  - `clawpack` (core claw, this is v4.6 (I think) of Clawpack)
  - `geoclaw`
      - The core of D-Claw is located in `D-Claw/geoclaw/2d/lib_dig`

- `python`, which contains two unpackaged python modules and some additional tools. These are made accessible to the runs by setting the `PYTHONPATH` to include them.

  - `dclaw`. Upon initial review, this is mostly code to deal with D-Claw specific pre or postprocessing. It has quite a few comments that are helpful to understand the D-Claw inputs and outputs.

  - `pyclaw` The python wrapping of clawpack. This provides the ability to specify and run clawpack solvers. The modern (clawpack v5.6) version has a packaged version of pyclaw, compatible with python 3.6.
