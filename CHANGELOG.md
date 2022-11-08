# Changelog

## known issues

1. `MKL` includes outdated `FEAST`, the external names in `FEAST` library are modified to avoid linking error.
2. `OpenBLAS` causes SEGFAULT with version 0.3.15+ when compiled with `DYNAMIC_ARCH` enabled.

## version 2.7

1. add `-nu` (`--noupdate`) flag to skip check of new version on startup
2. fix issue [#74](https://github.com/TLCFEM/suanPan/issues/74)

## version 2.6

1. update `MKL` to version `2022.2.0`
2. update `TBB` to version `2021.7.0`
3. update `VTK` to version `9.2.2`
4. add docker images and docker build scripts
5. add `TabularSpline` amplitude that uses cubic spline interpolation
6. add `upsampling` command to upsample time series data
7. add `sdof_response` command to compute response of single degree of freedom system
8. add `response_spectrum` command to compute response spectrum for given ground motion

## version 2.5

1. reformulate NM sections, add `NMB21E` element with end moment release
2. add couple stress membranes `CST3`, `CST6`, `CSM4-8`
3. add universal iterative solvers `BiCGSTAB` and `GMRES`, and preconditioners `Jacobi` and `ILU`
4. add support for `icx` and `ifx` compilers, add support for `clang` on linux
5. fix a bug in `GSSSS` with loads are applied as support motions, add `GSSSSOptimal` scheme
6. add `MassPoint2D` and `MassPoint3D` elements

## version 2.4

1. add `RestitutionWall` constraint which conserves momentum and energy
2. add `benchmark` command to benchmark platform
3. constraints and loads are processed in a fully parallelized manner
4. add 3D viscous damper `Damper03` and `Damper04`
5. bugfixes

## version 2.3

1. update `Armadillo` to version 11.0
2. relocate history record file under home folder
3. add `GSSSS` integrator
4. `LeeNewmark` now supports `PARDISO`, `CUDA` and `FGMRES` solvers
5. move to `C++20`, need `GCC 10.3.0`, `Clang 13.0.1`, `MSVC 14.31`
6. add `MOMENTUM` to record system momentum
7. use non-iterative algorithm for force based beams `F21`, `F21H` and `F31`

## version 2.2

1. add `example` command to showcase the creation of a simple model
2. update `VTK` to version 9.1.0
3. add `LogicAND`, `LogicOR` and `LogicXOR` convergers to use multiple criteria
4. update `oneMKL` to `2022.0.3` on Windows
5. move to `VS2022`

## version 2.1

1. update `Armadillo` to version 10.8
2. add recorder tag to recorded files, remove timestamp for hdf5 files, easier to manage different recorders
3. several minor bugfixes
4. improve `LeeNewmark` and `LeeNewmarkFull` performance
5. correct multithreaded `SuperLU` implementation, change default number of threads to 10
6. bugfixes regarding sparse matrix representation

## version 2.0

1. fix a bug in elastic stiffness in CDP model
2. add porous media plane strain elements `PCPE4UC`, `PCPE8UC`, `PCPE4DC`, `PCPE8DC`
3. add N-M interaction enabled beam element `NMB31` and `NMB21`
4. add N-M interaction enabled section `NM2D1`, `NM3D1` (elastic) and `NM2D2`, `NM3D2` (inelastic)
5. change Rayleigh damping related syntax to include tangent stiffness term
6. add different stiffness types to `LeeNewmarkFull`, add support of geometry nonlinearity
7. revise section definition
8. add `B3DC` corotational formulation support to 3D beams

## version 1.9

1. update `Armadillo` to version 10.7
2. update `TBB` version 2021.4.0
3. add `FGMRES` iterative solver
4. switch to `core20` on snap
5. fix the visualisation bug with installation via snap
6. add `LineUDL2D` and `LineUDL3D` loads

## version 1.8

1. add `PlaneSymmetric13` and `PlaneSymmetric23` wrappers
2. add `CoulombFriction` material
3. update `Armadillo` to version 10.6

## version 1.7

1. revise `SimpleSand` model
2. add `DafaliasManzari` sand model
3. add `materialtestbystrainhistory` and `materialtestbystresshistory` utility functions
4. bugfix: potential racing in initialising reference dof, change to serial initialisation
5. bugfix: wrong update of plastic strain in `CDP` model
6. add `CDPM2` model with isotropic damage

## version 1.6

1. add `terminal` command and `output_folder` setting
2. some minor updates

## version 1.5

1. add `scoop` support
2. add `Contact3D` 3D node-triangular facet contact element
3. add `NodeLine` and `NodeFacet` contact constraint
4. add `Sleeve2D`, `Sleeve3D`, `MaxGap2D` and `MaxGap3D` constraint
5. update `OpenBLAS` to version 0.3.15
6. update `Material` class to accommodate couple stress related quantities

## version 1.4

1. add `R2D2` and `R3D2` alias for fixed length constraint
2. add `MinGap2D` and `MinGap3D` inequality constraints
3. add `SupportMotion` loads, including `SupportDisplacement`, `SupportVelocity` and `SupportAcceleration`
4. update `Armadillo` to version 10.4
5. add functionality to check new version
6. add `FEAST` solver
7. improve mixed precision solver performance
8. add `CUDA` solver for dense matrix

## version 1.3

1. update handling of constraints and loads
2. store commands in backup file in CLI mode
3. add `FixedLength2D` and `FixedLength3D` nonlinear constraints
4. improve handling of constraints in dynamic analysis
5. add `NLE1D01` model
6. minor bugfixes

## version 1.2

1. remove dependency on `MAGMA`, now `CUDA` is directly used as the GPU solver
2. add `PARDISO` sparse solver
3. upgrade to `Intel oneAPI Toolkit`
4. add C interface material model
5. remove all reinforced elements, reinforcement can be handled by material models

## version 1.1

1. add phase field enabled elements: `DCP3`, `DCP4`, `DC3D4`, `DC3D8` elements
2. add support to record nodal damping/inertial force `DF` and `IF`
3. add regularized `Yeoh` model for compressible rubbers
4. improve stability of `RambergOsgood` model
5. add `LeeNewmarkFull` damping model, improve performance of `LeeNewmark` damping model
6. add shared memory `SuperLU` solver
7. add `Spike` solver for banded matrices
8. add displacement based beam element with end moment release: `B21EL` and `B21EH` elements
9. correct name of `Kelvin` model

## version 1.0

1. initial release
