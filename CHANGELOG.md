# Changelog

## known issues

1. Eigenanalysis does not support distributed computation.
2. Arc-length analysis is limited on both SMP and DMP, mainly due to the lack of determinant computation.

## version 3.8

1. add optional hourglassing control for `CP4R` and `C3D8R` elements
2. add `Prestrain` wrapper for uniaxial materials to apply prestrain [#266](https://github.com/TLCFEM/suanPan/pull/266)
3. add `GERKN` generalized explicit RKN time integration [#268](https://github.com/TLCFEM/suanPan/pull/268)
4. add `GSSE` explicit time integration [#274](https://github.com/TLCFEM/suanPan/pull/274)
5. fix improper detection of linear systems [#270](https://github.com/TLCFEM/suanPan/pull/270)
6. fix wrong corrector in explicit Bathe time integration [#272](https://github.com/TLCFEM/suanPan/pull/272)
7. improve numerical integration stability and robustness of various material models
8. update `OpenBLAS` to version `0.3.30`

## version 3.7

1. update `Armadillo` to version `14.4.1`
2. (breaking) remove general iterative solvers [#250](https://github.com/TLCFEM/suanPan/pull/250)
3. remove `MUMPS` and `lis` solvers for single-node binaries
4. update `MinGW-w64` with UCRT and GCC 13.3.0, see SDK [link](https://github.com/brechtsanders/winlibs_mingw/releases/download/13.3.0posix-11.0.1-ucrt-r1/winlibs-x86_64-posix-seh-gcc-13.3.0-mingw-w64ucrt-11.0.1-r1.7z)
5. add cluster support [#253](https://github.com/TLCFEM/suanPan/pull/253)
6. (breaking) refactor argument parser, some arguments are changed [#257](https://github.com/TLCFEM/suanPan/pull/257)

## version 3.6

1. add `AFCO1D` material with strain memory [#217](https://github.com/TLCFEM/suanPan/pull/217)
2. update `Catch2` to version `3.7.1`
3. add `Subloading1D` material [#219](https://github.com/TLCFEM/suanPan/pull/219)
4. add `Subloading~~Metal~~` material [#221](https://github.com/TLCFEM/suanPan/pull/221)
5. add `arm64` build
6. update `Armadillo` to version `14.2.3`
7. update `VTK` to version `9.4.0`
8. update `HDF5` to version `1.14.5`
9. update `OpenBLAS` to version `0.3.29`
10. revise US and EU section database
11. add support of `AMD Optimizing CPU Libraries (AOCL)` on linux
12. update `TBB` to version `2022.0.0`

## version 3.5

1. add `MaxForce` constraint [#204](https://github.com/TLCFEM/suanPan/pull/204)
2. update `Armadillo` to version `14.0.2`
3. update `OpenBLAS` to version `0.3.28`
4. update `Catch2` to version `3.7.0`
5. revise stiffness matrix formulation in shell elements [#208](https://github.com/TLCFEM/suanPan/pull/208)
6. add dev containers for easier DE setup

## version 3.4

1. update `Armadillo` to version `12.8.2` [#193](https://github.com/TLCFEM/suanPan/pull/193)
2. plane strain Duncan-Selig soil model [#195](https://github.com/TLCFEM/suanPan/pull/195)
3. update `Catch2` to version `3.5.4`
4. update `TBB` to version `2021.12.0` [#199](https://github.com/TLCFEM/suanPan/pull/199)
5. update `MUMPS` to version `5.7.1` [#202](https://github.com/TLCFEM/suanPan/pull/202)

## version 3.3

1. (breaking) revise syntax for `ConcreteTsai`, `Concrete21` and `Concrete22` using a more consistent definition
2. fix tangent stiffness in `ConcreteK4` model [#171](https://github.com/TLCFEM/suanPan/pull/171)
3. update `OpenBLAS` to version `0.3.25`
4. update `SuperLU` to version `6.0.1`
5. better `RCM` algorithm that may result in a smaller bandwidth, thus, potentially more efficient
   solving [#175](https://github.com/TLCFEM/suanPan/pull/175)
6. update `Armadillo` to version `12.6.7` [#180](https://github.com/TLCFEM/suanPan/pull/180)
7. enable lazy evaluation and avoid temporary global matrices, faster dynamic
   analysis [#183](https://github.com/TLCFEM/suanPan/pull/183)
8. bugfixes: [#185](https://github.com/TLCFEM/suanPan/pull/185)
9. update `Catch2` to version `3.5.2`
10. update `fmt` to version `10.2.1`
11. add nonviscous damping: `NonviscousNewmark` (global level integrator), `ElementalNonviscous` (element level
    modifier) and `Nonviscous01` (material level model)

## version 3.2

1. revise uniaxial universal damage models [#153](https://github.com/TLCFEM/suanPan/pull/153)
2. update `OpenBLAS` to version `0.3.24`
3. add a new uniaxial concrete model `ConcreteK4` [#155](https://github.com/TLCFEM/suanPan/pull/155)
4. add beam element for arbitrary thin-/thick-walled open/close section with torsion and
   warping `B31OS` [#159](https://github.com/TLCFEM/suanPan/pull/159)
5. better local iteration convergence criterion [#161](https://github.com/TLCFEM/suanPan/pull/161)
6. `B31OS` and `EB31OS` associated transformations `B3DOSL`, `B3DOSC`; sections, `Fibre3DOS`, `Cell3DOS`; material
   wrappers `OS146`, `OS146S`
7. add elemental damping using Lee's model
8. support Lode angle in CDPM2 [#163](https://github.com/TLCFEM/suanPan/pull/163)
9. add `AICN` cubic Newton solver [#165](https://github.com/TLCFEM/suanPan/pull/165)
10. remove `Bilinear2D` material, use `PlaneStress`/`PlaneStrain` wrapper and `BilinearJ2` 3D model instead

## version 3.1

1. iterative solvers by the Lis library [#145](https://github.com/TLCFEM/suanPan/pull/145)
2. update `Armadillo` to version `12.6.3` [#149](https://github.com/TLCFEM/suanPan/pull/149)
3. add `TimberPD` 3D material for timber [#151](https://github.com/TLCFEM/suanPan/pull/151)

## version 3.0

1. add experimental `MAGMA` based GPU sparse solver [#123](https://github.com/TLCFEM/suanPan/pull/123)
2. add nonlinear transformation for shell elements [#124](https://github.com/TLCFEM/suanPan/pull/124)
3. update `VTK` to version `9.2.6`
4. add `CustomNodeGroup` [#126](https://github.com/TLCFEM/suanPan/pull/126)
5. add `TranslationConnector` [#127](https://github.com/TLCFEM/suanPan/pull/127)
6. add `CustomAmplitude` [#129](https://github.com/TLCFEM/suanPan/pull/129)
7. update `Armadillo` to version `12.2` [#134](https://github.com/TLCFEM/suanPan/pull/134)
8. add `AsymmElastic1D` [#135](https://github.com/TLCFEM/suanPan/pull/135)
9. update `TBB` to version `2021.9.0`
10. update `MUMPS` to version `5.6.0`

## version 2.9

1. matrix optimisation
2. update `Catch2` to version `3.3.1`
3. update `TBB` to version `2021.8.0`
4. add mixed precision algorithm for `MUMPS` solver [#119](https://github.com/TLCFEM/suanPan/pull/119)
5. add `CustomDegradation`, `CustomGurson` and `CustomGurson1D` models
6. update `Armadillo` to version `12.0` [#121](https://github.com/TLCFEM/suanPan/pull/121)

## version 2.8

1. better on screen display with the `fmt` library [#99](https://github.com/TLCFEM/suanPan/pull/99)
2. add command `overview`
3. update `OpenBLAS` to version `0.3.21`
4. add Euler buckling load check for `T2D2` [#104](https://github.com/TLCFEM/suanPan/pull/104)
5. speed-up analysis with visualisation recorder [#102](https://github.com/TLCFEM/suanPan/pull/102)
6. update `VTK` to version `9.2.5`
7. add `Expression` to support custom function definition [#105](https://github.com/TLCFEM/suanPan/pull/105)
8. add `CustomMises1D`, `CustomCC`, `CustomCDP`, `CustomDP`, `CustomJ2` and `CustomHoffman` models

## version 2.7

1. optimise assembling of symmetric global matrices [#79](https://github.com/TLCFEM/suanPan/pull/79)
2. extend `BatheTwoStep` to allow customisation of spectral radius [#81](https://github.com/TLCFEM/suanPan/pull/81) and
   sub-step size [#82](https://github.com/TLCFEM/suanPan/pull/82)
3. update `Catch2` to version `2.13.10`
4. update `Armadillo` to version `11.4`
5. update modern `Arpack` [#94](https://github.com/TLCFEM/suanPan/pull/94)
6. add `Tchamwa` [#88](https://github.com/TLCFEM/suanPan/pull/88),
   `BatheExplicit` [#90](https://github.com/TLCFEM/suanPan/pull/90)
   and `GeneralisedAlphaExplicit` [#93](https://github.com/TLCFEM/suanPan/pull/93) explicit time integration methods
7. add `OALTS` two-step implicit time integration method [#92](https://github.com/TLCFEM/suanPan/pull/92)
8. add `Sinh1D` and `Tanh1D` nonlinear elastic 1D material
9. add `linear_system` flag to speed up linear system analysis

## version 2.6.1

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
5. (breaking) change Rayleigh damping related syntax to include tangent stiffness term
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
