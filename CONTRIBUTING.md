# Contribution

Pull requests are welcome.  
By submitting a pull request, you acknowledge and agree that project maintainers may modify your contribution without prior notice.
This may include changes to the submitted code and, where applicable, updates to the project's overall licensing terms.

If you intend to retain full copyright control over your work, please consider contributing through dynamically linked libraries.

Feature requests should be submitted by opening an issue.  
Templates are available for feature request submissions.

Private-domain collaboration in any form is also available; please contact `tlcfem(at)gmail(dot)com`.

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
