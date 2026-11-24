# Standalone FGMRES wrapper test

`test_fgmres.cpp` includes `Toolbox/fgmres.hpp` directly and links the real
oneMKL RCI FGMRES routines. It solves a small dense system, checks the residual
and iteration count, verifies that invalid input is rejected, and verifies that
an unattainable tolerance is reported as non-convergence instead of success.
No part of the suanPan application is compiled.

The `shims/` directory provides only the declarations of the four MKL entry
points and the two logging macros used by the wrapper; this avoids pulling the
full application header and its link dependencies into the test.

## Dependencies

A C++20 compiler and the oneMKL runtime are required. If oneMKL is not installed
system-wide, the Python wheel is enough for the runtime:

```bash
python3 -m venv /tmp/fgmres-venv
/tmp/fgmres-venv/bin/pip install mkl
```

## Build and run

```bash
MKL_LIB_DIR=/tmp/fgmres-venv/lib ./Toolbox/tests/fgmres/build.sh
```

`MKL_LIB_DIR` may be omitted when oneMKL is installed under
`/opt/intel/oneapi/mkl/latest/lib/intel64` or when running from a virtual
environment that has the `mkl` wheel installed.
