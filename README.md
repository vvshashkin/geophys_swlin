# Linearized shallow water model
## About
This model is intended to reproduce the following phenomenon from geophysical hydrodynamics
1. surface gravity waves,
1. surface inertia-gravity waves
1. Rossby-adjustment to the geostrophic ballance,
1. Rossby waves
1. equatorial wave dynamics
1. other phenomenon of linear shallow water dynamics on the rotating sphere.

The primary goal is demonstration for "Intriduction to geophysical hydrodynamics" lections.

## Installation and usage
### Compilation using Linux
To compile this model you basicly need two things:
1. fortran compiler (eg. gfortran),
1. make utility.
The second is available by default in any Linux distribution. The first can be installed using package manager. For instance in Ubuntu
```bash
    sudo apt install gfortran
```
works quite well.

To compile model clone the code for github and type
```bash
    make
```
in the code directory. This will create swlin executable

### Compilation using Windows
Recomendations for compiling model in windows are to follow.

### Running model
Once the model is compiled you can run it with typing
```bash
    ./swlin
```
in your command line in the model directory. Running the model will produce a number of 
```
 exp solver:          34  iterations, residual ~   5.2909890457487657E-013
 fields record          357 written
```
-like strings. The file f.dat will be created.

Model configuration is set using Fortran namelist. The model uses namelist.igw file by default. You can specify custom namelist in the forst command-line argument:
```bash
    ./swlin your_namelist_file_name
```

### Namelist guide
This is a typicall namelist for the model
```
$prm
NLON=240,
NLAT=120,
NSTEP = 360,
NZAP = 1,
DT = 900.,
$end

$cst
href=2000.
omega=2d-4,
lbeta=.true.,
$end

$ini
lam0 = 180.,
phi0 = 0.,
rhill = 1500d3,
hill_type="gauss",
$end
```
* NLON and NLAT are the number of points in longitude and latitude correspondingly
* NSTEP is the number of time steps to perform
* DT is the size of time-step in seconds
* NZAP is the number of time-steps passed between the solution writes to the disk (NSTEP=1 means the solution will be written to f.dat each time-step).

* href is the mean (background) fluid height
* omega is the planet rotation angular velocity
* lbeta turns on/off the latitudinal dependency of Corriolis parameter

* lam0 and phi0 are longitude and latitude position of initial height-field disturbance
* rhill is the height-field disturbance radius
* hill_type is the type of initial disturbance, possible options:
    * "gauss"=exp(-r^2/rhill^2),
    * "cone" - linear decay from center,
    * "cylinder" - 1 inside rhill, 0 outside,
    * "cosbell" - decay from center by trig. function


## Data visualization
Using [GRADS](http://cola.gmu.edu/grads/) is recomended. CTL-file f.ctl for grads is included in the model code.

## Model description
The model solves the following system of equations:

<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial\vec{v}}{\partial t}=-g\nabla h-f k \times \vec{v},">
<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial h}{\partial t}=-H_{0}\nabla \cdot \vec{v},">

where <img src="https://render.githubusercontent.com/render/math?math=\vec{v}"> is the horizontal wind vector, <img src="https://render.githubusercontent.com/render/math?math=h"> is fluid layer height perturbation, <img src="https://render.githubusercontent.com/render/math?math=H_0"> is the background (mean) fluid layer height.

We use 2-nd order accurate formulae for discretization of gradient and divergence operators on the regular latitude-longitude grid. C-staggered grid is used. Additional longitudinal filter (smoother) is applied to the fields in the polar latitudes to avoid instability. Matrix exponential method is used for time-integration. 
