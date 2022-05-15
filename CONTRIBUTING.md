# Contribution

Contributions are welcome. Feature requests are welcome, please email me at `tlcfem(at)gmail(dot)com`. If needed,
additional help/collaboration can be provided.

Code in high quality can be merged into the main repository. By doing so, the authors consent to modifications (
refactoring for maintenance, etc.) by other developers in future without notice.

## Code Format

Please use the file `Resharper.DotSettings` to format the code. It is the configuration that can be imported
into [ReSharper C++](https://www.jetbrains.com/resharper-cpp/).

The `.clang-format` generated can also be used to override format settings in various IDEs.

## Naming Convention

1. Please use `snake_case` instead of `CamelCase` for variables.
2. Please use `CamelCase` for class names. Abbreviations shall be in upper case, for example `NURBSSurface` instead
   of `NurbsSurface`.
3. Please use meaningful variable names instead of abstract names such as `a`, `tt`, etc.
4. Please avoid abbreviations, except for common ones such as `ptr` and `tmp`. For variables related to the
   corresponding theories, they can be spelled out. It does not hurt to define `epsilon_vol` instead of `ev` to
   represent volumetric strain variable. This is for the ease of readability and maintenance.

## Some Tips

1. Do not use variadic static variables.
2. Do not use raw pointers, instead, please use smart pointers.
3. For performance, do not use four (second) dimensional containers for fourth (second) order tensor if it is symmetric.
4. Use armadillo classes `mat` and `vec` in general cases. If not applicable, use STL containers `std::vector`.
5. Provide unit tests ([Catch2](https://github.com/catchorg/Catch2)) if applicable, otherwise provide an example model
   that tests the code and documentation.
6. Please check similar implemented code before getting started.
7. Apply early return.
8. Branch-free code is preferred.
