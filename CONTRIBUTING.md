# Contributing to suanPan

Thank you for your interest in **suanPan**. Contributions of all kinds are welcome: bug reports, feature requests, documentation, tests, example models and code.

This guide covers the development workflow, the coding style used in the project, and the release procedure.

## Contents

- [Before you contribute](#before-you-contribute)
- [Development workflow](#development-workflow)
- [Coding style](#coding-style)
- [Discouraged styles](#discouraged-styles)
- [Project layout](#project-layout)
- [Release procedure](#release-procedure)

## Before you contribute

1. Search the existing [issues](https://github.com/TLCFEM/suanPan/issues) and pull requests before starting.
2. Open an issue before large changes, breaking changes or significant refactors. This avoids duplicated work and lets the design be discussed first.
3. Keep each pull request focused. Do not mix unrelated refactoring, formatting changes and new features in one pull request.
4. Feature requests should be submitted by opening an issue. Templates are available for feature request submissions.

### Contribution terms

By submitting a pull request, you acknowledge and agree that project maintainers may modify your contribution without prior notice.
This may include changes to the submitted code and, where applicable, updates to the project's overall licensing terms.

If you intend to retain full copyright control over your work, please consider contributing through dynamically linked libraries.

Private-domain collaboration in any form is also available; please contact `tlcfem(at)gmail(dot)com`.

## Development workflow

The default development branch is `dev`. The `master` branch is used for releases.

1. Fork the repository and clone it with submodules:

   ```bash
   git clone --recursive https://github.com/TLCFEM/suanPan.git
   cd suanPan
   git switch dev
   ```

2. Create a topic branch from `dev`:

   ```bash
   git switch -c fix/my-change
   ```

3. Build and test your change. See the [manual](https://tlcfem.github.io/suanPan-manual/latest/Basic/Compile/) for the full dependency setup.

   ```bash
   cmake -S . -B build \
     -DCMAKE_BUILD_TYPE=Release \
     -DCMAKE_INSTALL_PREFIX=build/dist
   cmake --build build --target install --parallel
   ./build/dist/bin/suanPan -ct
   ```

   Unit tests are built into the main executable and are run with `-ct` (or `--catch2test`).

4. Format the files you changed:

   ```bash
   clang-format -i path/to/changed_file.cpp path/to/changed_file.h
   ```

   `./Script/Format.sh` formats the whole source tree. Use it only when you intend to reformat the whole tree, otherwise it creates unrelated diffs.

5. Commit with a short, imperative summary. Put the reason for the change and any important details in the commit body. Reference the relevant issue or pull request when applicable.

6. Open a pull request against `dev`, not `master`. Describe what changed, why it changed, how it was tested, and any breaking behaviour.

7. Make sure CI is green and address review comments. Maintainers may squash or edit commits when merging.

A pull request is ready to merge when it:

- builds without new warnings;
- passes the existing tests and includes tests for new behaviour;
- is formatted with the project `.clang-format`;
- contains no unrelated changes, commented-out code or debug output;
- documents user-visible changes, especially breaking changes.

## Coding style

### Formatting

Use the repository [`.clang-format`](.clang-format) file. It is based on the LLVM style with the following important settings:

- 4-space indentation, no tabs;
- no hard column limit;
- pointers and references attached to the type: `int* ptr`, `const std::string& name`;
- no space before parentheses: `if(...)`, `for(...)`, `function(...)`;
- opening braces on the same line for functions, classes and control statements; `else` and `catch` on a new line;
- includes sorted and regrouped;
- one space before trailing comments.

Also run [`.clang-tidy`](.clang-tidy) where practical. Do not disable checks or add blanket suppressions to silence warnings. Fix the cause, or suppress a single diagnostic with a short justification.

### Naming

| Item | Convention | Example |
| --- | --- | --- |
| Classes, structs, enums | `PascalCase`, acronyms fully capitalized | `NURBSSurface`, `CDPM2` |
| Functions, variables, members | `snake_case` | `compute_stiffness`, `elastic_modulus` |
| File-scope constants | `UPPER_SNAKE_CASE` | `MAX_ITERATION` |
| Test files | `Test<Subject>.cpp` | `TestMatrix.cpp` |
| Header/source files | Match the class or module name | `BilinearDP.h`, `BilinearDP.cpp` |

Additional naming rules:

1. Use meaningful names. Avoid abstract names such as `a`, `b`, `data2` when a descriptive name is possible.
2. Avoid abbreviations except common ones such as `ptr` and `tmp`.
3. Follow the naming already used by the surrounding module when extending existing code.

### Modern C++

The project uses C++20. Prefer modern, explicit and safe C++:

- Use RAII and smart pointers (`std::unique_ptr`, `std::make_unique`, `std::shared_ptr`, `std::make_shared`) for ownership.
- Use `nullptr`, `constexpr`, `enum class`, `using` and `static_cast`.
- Use `[[nodiscard]]` on getters and factory functions, and `override`/`final` on virtual functions.
- Use `const` correctness. Pass cheap types by value and non-trivial types by `const` reference.
- Use `auto` when the type is obvious from the initializer or when it avoids repetitive type names. Keep explicit types at public API boundaries.
- Use `std::array` and `std::vector` for general containers. Use Armadillo `mat`/`vec` (and the project aliases) for mathematical operations.
- Prefer early returns and guard clauses over deeply nested code.
- Branch-free code is preferred where it remains readable.
- Use the project logging macros (`suanpan_info`, `suanpan_debug`, `suanpan_warning`, `suanpan_error`, `suanpan_fatal`) instead of `std::cout` or `printf` in library code.
- Use the project's approximate comparison helpers for floating-point comparisons where appropriate.

### Constants and `constexpr`

Use `constexpr` for values known at compile time, and reserve plain `const` for runtime values and per-object instance state. Which specifiers to add depends on the scope:

| Scope | Preferred declaration |
| --- | --- |
| Static data member of a class or struct | `static constexpr`, initialised in the class; it is implicitly `inline` in C++20, so omit the out-of-line `.cpp` definition |
| Shared namespace-scope constant in a header | `inline constexpr`, so every translation unit refers to one entity with one address |
| Private namespace-scope constant in a `.cpp` file | plain `constexpr`; `static` is redundant because a namespace-scope `constexpr` already has internal linkage |
| Small scalar local in a non-`constexpr` function | plain `constexpr` for numbers, booleans and enums |
| Local array, lookup table or larger aggregate | `static constexpr`, so it is not materialised on the stack on every call |
| Local in a `constexpr` or `consteval` function | plain `constexpr` only; `static constexpr` is ill-formed there |
| `constexpr` reference | `inline constexpr` in a header, `constexpr` in a `.cpp` file, `static constexpr` in a class, always bound to an object with static storage duration |

Additional rules:

1. Use `constexpr` only when the type is a literal type and the initialiser is a constant expression. Armadillo `mat`, `vec` and `uvec`, as well as `TokenRegistry`, are not literal types; keep them as `static const` with an out-of-line definition (`inline static const` if they must be defined in a header).
2. Standard library functions such as `std::sqrt` and `std::log` are not `constexpr` in C++20, so constants initialised with them must remain `const`.
3. A non-`static` `const` data member such as `const double length = 0.;` is per-object instance state, often written through `access::rw`. Do not convert it into a `static constexpr` class member.
4. When a literal-type class member becomes `static constexpr`, delete its old out-of-line definition rather than keeping both declarations.

### Headers and includes

1. Keep the copyright and license header used by the file.
2. Use the existing include-guard style (`#ifndef ... #define ... #endif`). Do not introduce `#pragma once` in new files.
3. Include what you use. Do not rely on transitive includes and do not use `<bits/stdc++.h>`.
4. Do not put `using namespace ...;` in a header. Prefer fully qualified names or a local `using` declaration in a source file.
5. Keep the include order produced by `.clang-format`.
6. Add Doxygen comments for public classes and functions, following the existing `@class`, `@brief`, `@author`, `@date`, `@file` and `@addtogroup` style.

### Tests and examples

1. Provide [Catch2](https://github.com/catchorg/Catch2) tests for new behaviour when applicable. Add test files under `UnitTest/` and register them in `UnitTest/CMakeLists.txt`.
2. Tests must be deterministic. Avoid relying on random data unless the seed is fixed.
3. Use `REQUIRE`, `SECTION` and `Approx` in the same style as the existing tests.
4. When a unit test is not practical, provide example models under `Example/` or `Resource/` that cover the implementation as much as possible.
5. Check similar existing implementations before adding a new element, material, section, solver or command.

## Discouraged styles

| Avoid | Prefer |
| --- | --- |
| `CamelCase` variables or functions | `snake_case` |
| `NurbsSurface`, `HtmlParser` | `NURBSSurface`, `HTMLParser` |
| Owning raw pointers, `new` and `delete` | `std::unique_ptr`/`std::shared_ptr` with `std::make_unique`/`std::make_shared` |
| `NULL` or `0` as a pointer value | `nullptr` |
| C-style casts such as `(double)x` | `static_cast<double>(x)` |
| `typedef` | `using` |
| Mutable global state, non-const static variables or shared working buffers | Explicit parameters, local state or member state |
| `using namespace std;` in headers | Qualified names or a local `using` in a source file |
| Deeply nested `if`/`for` blocks | Early returns and guard clauses |
| Magic numbers | Named `constexpr` constants |
| `std::cout`, `printf` or `std::cerr` in library code | `suanpan_info`, `suanpan_warning`, `suanpan_error` |
| C-style arrays and manual bounds | `std::array`, `std::vector` or Armadillo types |
| Commented-out code and leftover debug prints | Remove them before opening the pull request |
| Catching everything with `catch(...)` and ignoring it | Catch the expected exception type; rethrow unknown exceptions |
| Unrelated reformatting or refactoring in a feature pull request | A separate, focused pull request |
| Abbreviations such as `stf`, `idx2`, `tmpmat` | Descriptive names such as `stiffness`, `second_index`, `temporary_matrix` |

## Project layout

- `CMakeLists.txt` — top-level CMake build file.
- `Include/` — public headers and vendored dependencies.
- `Libs/` — third-party or internal libraries.
- `UnitTest/` — Catch2 tests and the test harness.
- `Example/` — example models and usage samples.
- `Developer/` — templates for external elements, materials, modifiers and sections.
- `Script/` — helper scripts for build, format, versioning and packaging.
- `Resource/` — data and resource files used by examples and tests.
- `Domain/` — core domain objects such as nodes, elements and materials.
- `Solver/`, `Converger/`, `Constraint/`, `Element/`, `Material/`, `Section/`, `Checker/`, `Recorder/`, `Toolbox/`, `Step/`, `Load/`, `Database/` — core modules grouped by functionality.

When contributing, add or modify code in the appropriate module directory, update `CMakeLists.txt` as needed, and include tests or example models that demonstrate the behaviour. For large changes, open an issue first to discuss the design.

## Release procedure

The release procedure is for maintainers. It assumes write access to the code repository, the documentation repository, the Fedora rpkg/COPR project, the Snapcraft project and the Flathub repository.

Release tags in the code repository use the form `suanPan-vX.Y.Z`. Documentation tags use the version without the prefix, for example `4.3.1`.

### 1. Update the version and citation date

Run the version script from the repository root:

```bash
./Script/Version.sh 4 3 2
```

Then edit `CITATION.cff` and set `date-released` to the release date. Confirm that `version` in `CITATION.cff` matches the new version.

Add a new section at the top of `CHANGELOG.md` for the release. Describe user-visible changes, fixes and breaking changes, and link the relevant pull requests. Commit these changes on `dev`.

`Script/Version.sh` updates the version in the build files, packaging files and CI workflows. Review the diff to make sure every version reference is consistent.

### 2. Merge the release to `master`

Open a pull request from `dev` to `master` and merge it after CI passes. Merging to `master` triggers the `Stable` CI/CD workflow.

### 3. Promote the CI release

1. Open the Actions tab and wait for the `Stable` workflow for the release merge commit to finish. It creates a prerelease tagged `daily-YYMMDD` and uploads the platform archives.
2. Open the repository's Releases page and edit the `daily-YYMMDD` release.
3. Rename the tag to `suanPan-vX.Y.Z`.
4. Set the release title to `suanPan vX.Y.Z`.
5. Uncheck "Set as a pre-release" and check "Set as the latest release".
6. Replace the auto-generated body with the matching `CHANGELOG.md` section. Keep the VirusTotal results if they are useful.
7. Publish the release.

Do not create a separate release from scratch. Rename and edit the CI-created release so that the uploaded assets remain attached.

### 4. Trigger the Fedora rpkg build

1. Confirm that `Version:` in `suanPan.spec.rpkg` matches the new version.
2. Push the updated spec to the rpkg repository used by COPR.
3. Trigger a new build in the COPR project `tlcfem/suanPan`.
4. Verify that the build succeeds for all supported Fedora and EPEL chroots.

The spec uses the `rpkg` `git_dir_*` macros, so build from the release tag or the corresponding release commit.

### 5. Trigger the Snapcraft build

1. Confirm that `version` in `snapcraft.yaml` matches the new version.
2. Trigger new builds for `amd64` and `arm64` in the Snap Store.
3. After the builds pass review, promote the new revision to the stable channel.
4. Verify the published version with `snap info suanpan`.

### 6. Create the Flatpak release pull request

1. Run the `Flatpak` workflow from the Actions tab, or open the pull request manually.
2. The pull request must target `flathub/io.github.tlcfem.suanPan` on the `beta` branch.
3. Update the pinned source commit and the release date in `io.github.tlcfem.suanPan.yml`.
4. Confirm that the pinned commit is the release commit, not an arbitrary `dev` commit.
5. Wait for Flathub CI to pass and merge the pull request.

### 7. Release the documentation

1. Open the documentation repository, `TLCFEM/suanPan-manual`.
2. Merge the documentation changes into `master`.
3. Confirm that Read the Docs rebuilds the `latest` version successfully.

### 8. Tag the documentation repository

Create a tag without the `suanPan-v` prefix and push it:

```bash
git switch master
git pull
git tag 4.3.2
git push origin 4.3.2
```

The tag triggers the Read the Docs build for that version. Confirm that the new version appears in the Read the Docs version selector and that the build succeeds.

### Release checklist

- [ ] Version updated with `Script/Version.sh`
- [ ] `CITATION.cff` date and version updated
- [ ] `CHANGELOG.md` entry added
- [ ] Release merged to `master`
- [ ] `Stable` workflow completed successfully
- [ ] `daily-YYMMDD` release renamed to `suanPan-vX.Y.Z` with release notes
- [ ] Fedora rpkg/COPR build triggered and green
- [ ] Snapcraft builds triggered and promoted
- [ ] Flatpak pull request merged
- [ ] Documentation merged to `master`
- [ ] Documentation tag pushed and Read the Docs build green
