# Changelog

## version 1.x

1. update handling of constraints
2. store commands in backup file in CLI mode
3. add `FixedLength2D` and `FixedLength3D` nonlinear constraints

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
