#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# include_tree.py - Build and clean up the C/C++ include graph
#
# The project contains a number of first party source files and a few bundled
# third party libraries (Include/ and Toolbox/*-src/).  Only first party files
# are considered here; the bundled libraries are treated as external leaves.
#
# The script can
#   1. print the include tree of one or more files,
#   2. report duplicated/redundant includes,
#   3. remove only the duplicates that are guaranteed to be safe on every
#      platform and in every configuration,
#   4. replace standard C headers with their C++ counterparts.
#
# Include directives are parsed together with their preprocessor context.  An
# include is only considered redundant when the header that provides it is
# reachable through an *unconditional* chain of includes.  This is important:
# for example Toolbox/container.h provides <set> only in its non-TBB branch, so
# a file that uses std::set cannot rely on container.h in a SUANPAN_MT build.
#
# Usage:
#   python3 Script/include_tree.py --report
#   python3 Script/include_tree.py --tree Element/Beam/EB21.cpp
#   python3 Script/include_tree.py --all-trees --output /tmp/include_tree.txt
#   python3 Script/include_tree.py --fix --fix-stl
#   python3 Script/include_tree.py --fix --format
# -----------------------------------------------------------------------------

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Iterator

ROOT = Path(__file__).resolve().parent.parent

HEADER_EXTENSIONS = (".h", ".hh", ".hpp", ".hxx")
SOURCE_EXTENSIONS = (".c", ".cc", ".cpp", ".cxx")
ALL_EXTENSIONS = HEADER_EXTENSIONS + SOURCE_EXTENSIONS

# Directories that contain dependencies or build output rather than project
# source files.  Toolbox/*-src is handled by the glob below.
EXCLUDED_TOP_LEVEL = {"Include", "Libs", "build", "cmake-build-debug"}
EXCLUDED_ANYWHERE = {"__pycache__"}

# Standard C++ headers.  These are platform agnostic and therefore safe to
# prune when another kept include provides them unconditionally.
STL_HEADERS = frozenset(
    {
        "algorithm", "any", "array", "atomic", "barrier", "bit", "bitset",
        "cassert", "cctype", "cerrno", "cfenv", "cfloat", "charconv", "chrono",
        "cinttypes", "climits", "clocale", "cmath", "codecvt", "compare",
        "complex", "concepts", "condition_variable", "coroutine", "csetjmp",
        "csignal", "cstdarg", "cstddef", "cstdint", "cstdio", "cstdlib",
        "cstring", "ctime", "cuchar", "cwchar", "cwctype", "deque", "exception",
        "execution", "expected", "filesystem", "flat_map", "flat_set", "format",
        "forward_list", "fstream", "functional", "future", "generator",
        "initializer_list", "iomanip", "ios", "iosfwd", "iostream", "istream",
        "iterator", "latch", "limits", "list", "locale", "map", "mdspan",
        "memory", "memory_resource", "mutex", "new", "numbers", "numeric",
        "optional", "ostream", "print", "queue", "random", "ranges", "ratio",
        "regex", "scoped_allocator", "semaphore", "set", "shared_mutex",
        "source_location", "span", "sstream", "stack", "stdexcept", "stdfloat",
        "stop_token", "streambuf", "string", "string_view", "strstream",
        "syncstream", "system_error", "thread", "tuple", "type_traits",
        "typeindex", "typeinfo", "unordered_map", "unordered_set", "utility",
        "valarray", "variant", "vector", "version",
    }
)

# Prefer modern C++ headers over the C compatibility headers.
C_HEADER_MAP = {
    "assert.h": "cassert",
    "complex.h": "complex",
    "ctype.h": "cctype",
    "errno.h": "cerrno",
    "fenv.h": "cfenv",
    "float.h": "cfloat",
    "inttypes.h": "cinttypes",
    "limits.h": "climits",
    "locale.h": "clocale",
    "math.h": "cmath",
    "setjmp.h": "csetjmp",
    "signal.h": "csignal",
    "stdarg.h": "cstdarg",
    "stddef.h": "cstddef",
    "stdint.h": "cstdint",
    "stdio.h": "cstdio",
    "stdlib.h": "cstdlib",
    "string.h": "cstring",
    "time.h": "ctime",
    "uchar.h": "cuchar",
    "wchar.h": "cwchar",
    "wctype.h": "cwctype",
}

INCLUDE_RE = re.compile(r'^\s*#\s*include\s*([<"])([^>"]+)[>"]')
DIRECTIVE_RE = re.compile(r"^\s*#\s*(if|ifdef|ifndef|elif|else|endif)\b\s*(.*)")
DEFINE_RE = re.compile(r"^\s*#\s*define\s+([A-Za-z_]\w*)\b")
NOT_DEFINED_RE = re.compile(r"^!\s*defined\s*\(?\s*([A-Za-z_]\w*)\s*\)?$")

Condition = tuple[tuple[str, str], ...]


@dataclass(frozen=True)
class Include:
    name: str
    start: int
    end: int
    line: int
    condition: Condition = ()
    angled: bool = True

    @property
    def key(self) -> str:
        return self.name


@dataclass
class Duplicate:
    path: Path
    include: Include
    kind: str
    provider: str
    safe: bool
    reason: str = ""

    def describe(self) -> str:
        line = f"{self.path.relative_to(ROOT)}:{self.include.line}"
        return f"[{self.kind}] {line}: {self.include.name} <- {self.provider} ({self.reason})"


def is_excluded(relative: Path) -> bool:
    parts = relative.parts
    if not parts:
        return True
    if parts[0] in EXCLUDED_TOP_LEVEL or parts[0].startswith("build"):
        return True
    if parts[0].startswith("cmake-build") or parts[0].startswith(".git"):
        return True
    if any(part.startswith(".") for part in parts[:-1]):
        return True
    if any(part in EXCLUDED_ANYWHERE for part in parts[:-1]):
        return True
    # Third party source trees such as Toolbox/superlu-src.
    if any(part.endswith("-src") for part in parts[:-1]):
        return True
    return False


def collect_sources(root: Path = ROOT) -> list[Path]:
    files: list[Path] = []
    for path in root.rglob("*"):
        if not path.is_file() or path.suffix not in ALL_EXTENSIONS:
            continue
        relative = path.relative_to(root)
        if is_excluded(relative):
            continue
        files.append(path.resolve())
    return sorted(files)


def _is_guard(lines: list[str], index: int, rest: str) -> bool:
    """Return True if the directive at index is the start of an include guard."""
    match = DIRECTIVE_RE.match(lines[index])
    if not match:
        return False
    directive = match.group(1)
    macro: str | None = None
    if directive == "ifndef":
        parts = rest.split()
        if parts:
            macro = parts[0]
    elif directive == "if":
        not_defined = NOT_DEFINED_RE.match(rest)
        if not_defined:
            macro = not_defined.group(1)
    if not macro:
        return False
    for line in lines[index + 1 :]:
        stripped = line.strip()
        if not stripped or stripped.startswith("//") or stripped.startswith("/*") or stripped.startswith("*"):
            continue
        define = DEFINE_RE.match(line)
        return bool(define and define.group(1) == macro)
    return False


def parse_includes(path: Path) -> list[Include]:
    """Parse includes and record the active non-guard preprocessor condition."""
    lines = path.read_text(encoding="utf-8", errors="surrogateescape").splitlines()
    includes: list[Include] = []
    stack: list[tuple[str, str]] = []
    index = 0
    while index < len(lines):
        line = lines[index]
        directive = DIRECTIVE_RE.match(line)
        if directive:
            kind, rest = directive.group(1), directive.group(2).strip()
            if kind in ("if", "ifdef", "ifndef"):
                # Include guards (classic #ifndef X / #define X or #if !defined)
                # are always active, so they must not make an include
                # conditional.
                if not stack and _is_guard(lines, index, rest):
                    stack.append(("guard", rest))
                else:
                    stack.append((kind, rest))
            elif kind == "elif":
                if stack:
                    stack[-1] = ("elif", rest)
            elif kind == "else":
                if stack:
                    stack[-1] = ("else", "")
            elif kind == "endif":
                if stack:
                    stack.pop()
            index += 1
            continue

        match = INCLUDE_RE.match(line)
        if match:
            condition = tuple(item for item in stack if item[0] != "guard")
            includes.append(
                Include(
                    name=match.group(2),
                    start=index,
                    end=index + 1,
                    line=index + 1,
                    condition=condition,
                    angled=match.group(1) == "<",
                )
            )
        index += 1
    return includes


class IncludeGraph:
    def __init__(self, root: Path = ROOT) -> None:
        self.root = root
        self.files = collect_sources(root)
        self.valid = set(self.files)
        self.includes = {path: parse_includes(path) for path in self.files}
        self.include_dirs = [
            root,
            root / "Include",
            root / "Include" / "fmt" / "include",
        ]

    def resolve(self, include: Include, current: Path) -> Path | None:
        candidates: list[Path] = []
        if not include.angled:
            candidates.append(current.parent / include.name)
        candidates.extend(directory / include.name for directory in self.include_dirs)
        for candidate in candidates:
            try:
                candidate = candidate.resolve()
            except OSError:
                continue
            if candidate.is_file():
                return candidate
        return None

    def is_project(self, path: Path | None) -> bool:
        return path is not None and path in self.valid

    def project_children(self, path: Path) -> Iterator[tuple[Include, Path]]:
        for include in self.includes.get(path, []):
            child = self.resolve(include, path)
            if self.is_project(child):
                yield include, child  # type: ignore[arg-type]

    def external_children(self, path: Path) -> Iterator[Include]:
        for include in self.includes.get(path, []):
            if not self.is_project(self.resolve(include, path)):
                yield include

    @lru_cache(maxsize=None)
    def unconditional_closure(self, path: Path) -> frozenset[Path]:
        """Project headers reachable through unconditional includes only."""
        seen: set[Path] = set()
        stack = [path]
        while stack:
            current = stack.pop()
            for include, child in self.project_children(current):
                if include.condition:
                    continue
                if child not in seen:
                    seen.add(child)
                    stack.append(child)
        return frozenset(seen)

    @lru_cache(maxsize=None)
    def unconditional_external_names(self, path: Path) -> frozenset[str]:
        """External (system/third party) headers reachable without conditions."""
        names: set[str] = set()
        seen: set[Path] = set()
        stack = [path]
        while stack:
            current = stack.pop()
            if current in seen:
                continue
            seen.add(current)
            for include, child in self.project_children(current):
                if include.condition:
                    continue
                stack.append(child)
            for include in self.external_children(current):
                if include.condition:
                    continue
                names.add(include.name)
        return frozenset(names)

    def companion_header(self, source: Path) -> Path | None:
        for extension in (".h", ".hpp", ".hh", ".hxx"):
            candidate = source.with_suffix(extension)
            if candidate in self.valid:
                return candidate
        return None

    def find_duplicates(self, stl_in_headers: bool = False) -> list[Duplicate]:
        duplicates: list[Duplicate] = []
        for path in self.files:
            includes = self.includes[path]

            # 1. The exact same header included twice in the same condition.
            #    Includes under different conditions are kept because they may
            #    be needed by different platforms/configurations.
            groups: dict[tuple[str, Condition], list[Include]] = defaultdict(list)
            for include in includes:
                resolved = self.resolve(include, path)
                key = str(resolved) if resolved else include.name
                groups[(key, include.condition)].append(include)
            for (_, condition), group in groups.items():
                if len(group) < 2:
                    continue
                for include in group[1:]:
                    duplicates.append(
                        Duplicate(
                            path,
                            include,
                            "exact",
                            str(path.relative_to(ROOT)),
                            True,
                            "repeated directive" if not condition else "repeated directive in the same branch",
                        )
                    )

            # 2. An implementation file repeating a header that is already
            #    provided by its own (first) header.  The own header is always
            #    included first and the path is unconditional, so the repeated
            #    include can never be needed on any platform or configuration.
            header = self.companion_header(path)
            if header is not None and path.suffix in SOURCE_EXTENSIONS:
                own = [
                    include
                    for include in includes
                    if self.resolve(include, path) == header and not include.condition
                ]
                if own:
                    provided = self.unconditional_closure(header)
                    for include in includes:
                        resolved = self.resolve(include, path)
                        if resolved is None or resolved == header or resolved not in provided:
                            continue
                        duplicates.append(
                            Duplicate(
                                path,
                                include,
                                "own-header",
                                str(header.relative_to(ROOT)),
                                True,
                                "already included by the implementation's own header",
                            )
                        )

            # 3. A project header or standard header that another kept include
            #    already provides unconditionally.  This is the same idea as
            #    the own-header rule above, but uses any sibling include.  The
            #    provider must include the target without preprocessor
            #    conditions, otherwise the target may be needed in a
            #    configuration in which the provider does not expose it.
            providers: list[tuple[Include, frozenset[Path], frozenset[str]]] = []
            for include in includes:
                if include.condition:
                    continue
                resolved = self.resolve(include, path)
                if self.is_project(resolved):
                    providers.append(
                        (
                            include,
                            self.unconditional_closure(resolved),  # type: ignore[arg-type]
                            self.unconditional_external_names(resolved),  # type: ignore[arg-type]
                        )
                    )
            for include in includes:
                if include.condition:
                    continue
                resolved = self.resolve(include, path)
                is_project = self.is_project(resolved)
                if not is_project and (path.suffix not in SOURCE_EXTENSIONS and not stl_in_headers):
                    continue
                if not is_project and include.name not in STL_HEADERS:
                    continue
                for provider, project_names, external_names in providers:
                    if provider.line == include.line:
                        continue
                    if is_project:
                        if resolved in project_names:
                            duplicates.append(
                                Duplicate(
                                    path,
                                    include,
                                    "cross-project",
                                    provider.name,
                                    True,
                                    "already included unconditionally by another include",
                                )
                            )
                            break
                    elif include.name in external_names:
                        duplicates.append(
                            Duplicate(
                                path,
                                include,
                                "stl",
                                provider.name,
                                True,
                                "standard header provided unconditionally by another include",
                            )
                        )
                        break

        # De-duplicate identical findings.  A header may match several rules;
        # report it once, preferring the most specific/strongest rule.
        priority = {"exact": 0, "own-header": 1, "cross-project": 2, "stl": 3}
        unique: dict[tuple[Path, int], Duplicate] = {}
        for duplicate in duplicates:
            key = (duplicate.path, duplicate.include.line)
            previous = unique.get(key)
            if previous is None or priority.get(duplicate.kind, 99) < priority.get(previous.kind, 99):
                unique[key] = duplicate
        return sorted(unique.values(), key=lambda item: (str(item.path), item.include.line, item.kind))

    def modernize_c_headers(self) -> list[tuple[Path, Include, str]]:
        """Find standard C headers that should be replaced by C++ headers."""
        findings: list[tuple[Path, Include, str]] = []
        for path in self.files:
            for include in self.includes[path]:
                replacement = C_HEADER_MAP.get(include.name)
                if replacement is not None:
                    findings.append((path, include, replacement))
        return findings

    def tree(self, path: Path) -> str:
        lines: list[str] = []

        def walk(current: Path, prefix: str, ancestors: set[Path]) -> None:
            includes = self.includes.get(current, [])
            for index, include in enumerate(includes):
                last = index == len(includes) - 1
                branch = "└── " if last else "├── "
                condition = ""
                if include.condition:
                    condition = "  [" + " && ".join(expr for _, expr in include.condition) + "]"
                child = self.resolve(include, current)
                if self.is_project(child):
                    marker = " (already shown)" if child in ancestors else ""  # type: ignore[operator]
                    lines.append(f"{prefix}{branch}{include.name}{condition}{marker}")
                    if child not in ancestors:
                        walk(child, prefix + ("    " if last else "│   "), ancestors | {child})  # type: ignore[operator]
                else:
                    lines.append(f"{prefix}{branch}{include.name}{condition} [external]")

        lines.append(str(path.relative_to(ROOT)))
        walk(path, "", {path})
        return "\n".join(lines)


def _remove_lines(path: Path, include: Include) -> None:
    text = path.read_text(encoding="utf-8", errors="surrogateescape")
    lines = text.splitlines(keepends=True)
    del lines[include.start : include.end]
    path.write_text("".join(lines), encoding="utf-8", errors="surrogateescape")


def _replace_include(path: Path, include: Include, replacement: str) -> None:
    text = path.read_text(encoding="utf-8", errors="surrogateescape")
    lines = text.splitlines(keepends=True)
    line = lines[include.start]
    match = INCLUDE_RE.match(line)
    if not match:
        return
    opening, closing = match.group(1), ">" if match.group(1) == "<" else '"'
    lines[include.start] = INCLUDE_RE.sub(
        f'#include {opening}{replacement}{closing}',
        line,
        count=1,
    )
    path.write_text("".join(lines), encoding="utf-8", errors="surrogateescape")


def describe_duplicates(duplicates: Iterable[Duplicate], verbose: bool = False) -> str:
    counts: dict[str, int] = defaultdict(int)
    for duplicate in duplicates:
        counts[duplicate.kind] += 1
    lines = [f"Duplicate include findings: {len(duplicates)}"]
    for kind in sorted(counts):
        lines.append(f"  {kind}: {counts[kind]}")
    if duplicates:
        lines.append("")
        for duplicate in duplicates:
            if not verbose and not duplicate.safe:
                continue
            lines.append("  " + duplicate.describe())
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=ROOT, help="repository root")
    parser.add_argument("--tree", action="append", default=[], metavar="FILE", help="print the include tree of FILE")
    parser.add_argument("--all-trees", action="store_true", help="print the include tree of every scanned file")
    parser.add_argument("--report", action="store_true", help="report duplicated includes (default)")
    parser.add_argument("--fix", action="store_true", help="remove safely duplicated includes")
    parser.add_argument("--fix-stl", action="store_true", help="also remove duplicated standard headers from source files")
    parser.add_argument("--stl-in-headers", action="store_true", help="allow --fix-stl to modify headers as well")
    parser.add_argument("--modernize-c-headers", action="store_true", help="replace standard C headers with C++ headers")
    parser.add_argument("--format", action="store_true", help="run Script/Format.sh after changing files")
    parser.add_argument("--output", type=Path, help="write the report/trees to a file instead of stdout")
    parser.add_argument("--verbose", action="store_true", help="show non-removable findings as well")
    args = parser.parse_args(argv)

    root = args.root.resolve()
    graph = IncludeGraph(root)

    output: list[str] = []
    target_files: list[Path] = []
    for requested in args.tree:
        candidate = (root / requested).resolve()
        if candidate not in graph.valid:
            matches = [path for path in graph.files if path.name == requested or str(path).endswith(requested)]
            if len(matches) != 1:
                parser.error(f"cannot identify tree target: {requested}")
            candidate = matches[0]
        target_files.append(candidate)
    if args.all_trees:
        target_files = graph.files

    for target in target_files:
        output.append(graph.tree(target))
        output.append("")

    duplicates = graph.find_duplicates(stl_in_headers=args.stl_in_headers)
    if args.report or not (args.fix or args.tree or args.all_trees):
        output.append(describe_duplicates(duplicates, verbose=args.verbose))

    c_headers = graph.modernize_c_headers()
    if c_headers:
        output.append("")
        output.append(f"Standard C headers that can be replaced by C++ headers: {len(c_headers)}")
        for path, include, replacement in c_headers:
            output.append(f"  {path.relative_to(root)}:{include.line}: {include.name} -> {replacement}")

    def removable_items(items: list[Duplicate]) -> list[Duplicate]:
        selected = [item for item in items if item.safe and (item.kind != "stl" or args.fix_stl)]
        # The same directive may be reported by more than one rule.
        return list({(item.path, item.include.line): item for item in selected}.values())

    if args.fix:
        total_removed = 0
        total_changed: set[Path] = set()
        passes = 0
        # Removing one include can make another include redundant, so repeat
        # until the graph reaches a fixed point.  A hard limit protects
        # against a malformed graph.
        while True:
            removable = removable_items(duplicates)
            replacements = c_headers if args.modernize_c_headers else []
            if not removable and not replacements:
                break
            if passes >= 20:
                output.append("Stopped after 20 passes without reaching a fixed point.")
                break
            passes += 1
            changed: set[Path] = set()
            # Remove from the bottom up so earlier line numbers stay valid.
            for item in sorted(removable, key=lambda entry: (str(entry.path), -entry.include.line)):
                _remove_lines(item.path, item.include)
                changed.add(item.path)
            for path, include, replacement in replacements:
                _replace_include(path, include, replacement)
                changed.add(path)
            total_removed += len(removable)
            total_changed |= changed
            output.append("")
            output.append(f"Pass {passes}: removed {len(removable)} duplicated include(s) from {len(changed)} file(s).")
            if not changed:
                break
            graph = IncludeGraph(root)
            duplicates = graph.find_duplicates(stl_in_headers=args.stl_in_headers)
            c_headers = graph.modernize_c_headers()

        output.append("")
        output.append(f"Removed {total_removed} duplicated include(s) from {len(total_changed)} file(s) in {passes} pass(es).")
        if total_changed:
            output.append("")
            output.extend(f"  {path.relative_to(root)}" for path in sorted(total_changed))
        if total_changed and args.format:
            formatter = Path(root) / "Script" / "Format.sh"
            subprocess.run([str(formatter)], cwd=root, check=True)
            output.append("")
            output.append("Formatted with Script/Format.sh.")
        final_duplicates = graph.find_duplicates(stl_in_headers=args.stl_in_headers)
        output.append("")
        output.append(f"Remaining removable includes after multipass: {len(removable_items(final_duplicates))}.")
    else:
        removable = removable_items(duplicates)
        if removable:
            output.append("")
            output.append(f"Re-run with --fix to remove {len(removable)} safely duplicated include(s).")

    report = "\n".join(output).rstrip() + "\n"
    if args.output:
        args.output.write_text(report, encoding="utf-8")
    else:
        sys.stdout.write(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
