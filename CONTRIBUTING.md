# Contribution

Pull requests are welcome.
By submitting pull requests, consent to future modifications without notice is implied.

Feature requests can be made by creating issues.
There is a template for feature requests.

Collaborations of any form in the private domain are possible, please contact me at `tlcfem(at)gmail(dot)com`.

## Code Format Style

Please use the file `.clang-format` to format the code.

## Naming Convention

1. Use `snake_case` instead of `CamelCase` for variables.
2. Use `CamelCase` for class names. Acronyms shall be capitalized, for example `NURBSSurface` instead
   of `NurbsSurface`.
3. Use meaningful variable names instead of abstract names if possible.
4. Avoid abbreviations, except for common ones such as `ptr` and `tmp`, for the ease of readability and maintenance.

## Miscellaneous Tips

1. Avoid non-const static variables, avoid shared working buffers.
2. Use smart pointers instead of raw pointers.
3. Use armadillo classes `mat` and `vec` for variables that involve mathematical operations.
4. Use STL containers `std::array` and `std::vector` in general cases.
5. Provide [Catch2](https://github.com/catchorg/Catch2) tests if applicable, otherwise provide example models
   that cover the implementation as much as possible.
6. Check similar existing implemented code first.
7. Apply early return if possible, avoid deep nesting of code.
8. Branch-free code is preferred.
