# Style Guide
A variety of signature features of C++ have been used in Tinker9 code base:
- procedural programming, C and improved C features;
- object-oriented programming;
- generic programming (template);
- the C++ standard template library,

which in my opinion are four different programming languages. The language’s
great flexibility is intrinsically coupled with great complexity that can be
greatly eased by following a brief coding style.

## Table of Contents
- [Source Files and Format](#format)
- [C++ Syntax Version](#cpp.ver)
- [Global Variables](#glob.var)
- [Naming](#name)
   - [Constant Names](#name.const)
   - [Macro Names](#name.macro)
   - [Type Names](#name.type)
   - [Function Names](#name.func)

<a name='format'></a>
## Source Files and Format
Source code must be formatted before committed into the Git history. The file
types in the table below should be processed by
[clang-format.](https://clang.llvm.org/docs/ClangFormat.html)

| File Types | Filename Extensions |
|------------|---------------------|
| C/C++/CUDA header                    | `.h`, `.hh`           |
| C source                             | `.c`                  |
| C++ source (generic)                 | `.cpp`, `.cxx`, `.cc` |
| C++ source (with OpenACC directives) | `.cpp`                |
| CUDA source                          | `.cu`                 |

> Generic C++ files should not have any OpenACC directive or CUDA specific macro
> or syntax, such as `__global__` and triple arrow brackets `<<<>>>`.

The minimum version of clang-format required is 10. The following command will
read the configurations from `tinker9/_clang-format` (enabled by `-style=file`)
then format the source code in-place (`-i`).

```bash
# to format "include/ebond.h"
clang-format -i -style=file include/ebond.h
```

One limitation of clang-format is it doesn’t format `#pragma` properly at this
point, especially for the long directives. The way I am handling this issue now,
is first to replace `#pragma` to `//#prag`, then format the source code, and
change it back to `#pragma` in the end.

```bash
# my custom code-format.sh
# on Linux set SED to "sed -i"
# on macOS set SED to "perl -i -pe"

SED 's/#pragma /\/\/#prag /g' “$*”
clang-format -i -style=file “$*”
SED 's/\/\/ *#prag /#pragma /g' “$*”
```

We are open to suggestions on ways to format Python code.

<a name='cpp.ver'></a>
## C++ Syntax Version
C++ syntax should target C++11 standard.

C++14 and C++17 are good, but the support from the toolchain to the newer
standards may be years behind. As we are using multiple compilers, the
intersection of the commonly supported syntax features may not leave us many
options. And more likely, the options are limited by the combination of several
other problems. For instance, the NVIDIA driver is old and cannot be updated
because the GPU is old and you don’t have root privilege of the cluster nodes,
so CUDA is limited to 10.0, which only supports C++14 but not C++17. There are
also a lot of machines still using GCC 4.8 as of 2020, which doesn’t have a
complete C++11 support.

<a name='glob.var'></a>
## Global Variables
Depending on different situations, we use different ways to declare new global
variables.

If you want to provide a new global variable as a public interface similar to a
new variable inside a Fortran module, you should follow the standard C/C++
method: declaring it in a header file and define it in a source file.

```cpp
// new_var.h
extern int new_var_unstable_name;

// new_var.cpp
int new_var_unstable_name;
```

If the name and usage have been stable for a while, you should change it to the
following style.

```cpp
// new_module.h
TINKER_EXTERN int new_var_stable;
```

On the other hand, if you are adding an internal global variable, for instance a
preallocated global work array, you should always declare it inside the `detail`
namespace.

```cpp
// new_array.h
namespace detail {
extern int* global_work_array;
}

// new_array.cpp
namespace detail {
int* global_work_array;
}
```

<a name='name'></a>
## Naming
We always put code readability first. Abbreviations are acceptable if the
readability is not compromised, in situations where they are:
- widely used, e.g., `PBC`;
- straightforward with no ambiguity, e.g., `config` for configuration;
- names of methods in the literature and documented in the source code, e.g.,
  `GK` for AMOEBA Generalized Kirkwood Model.

There are three common styles we are using: `snake_case` `UpperCamelCase` (aka
`PascalCase`), and `ALL_CAPS`.

<a name='name.const'></a>
### Constant Names
Constants in the source code usually appear as `enum`, `class enum`, and
`constexpr` values. Avoid defining constants in the macros unless you have no
alternative options. All caps is usually the preferred style, but snake case is
also acceptable (for instance, set `constexpr float sqrt_2pi = 2.50662827463;`
in the math calculation).

<a name='name.macro'></a>
### Macro Names
Although there are a few exceptions in the code base for the function-like
macros, most macro names are all caps. Always choose constant variables over
macros.

<a name='name.type'></a>
### Type Names
Light-weight types (e.g., `real`) prefer snake case. Heavy-weighted types, such
as complicated class names, should always use upper camel case.

<a name='name.func'></a>
### Function Names
We use snake case for function names in general and we follow a set of
conventions to name functions.

By default, we assume every function is a regular CPU function unless it has
been explicitly marked. For example, in the `energy()` function where we call
`elj()` to calculate Lennard-Jones energy, we should not care about or have any
assumption on which hardware (CPU vs. GPU) or code path (OpenACC vs. CUDA) would
be used.

Inside function `elj()`, we would then have to choose among different possible
code paths. We use **different suffix** to indicate different **CPU function**
that has access to different internal implementation:
- `elj_acc()`: OpenACC implementation;
- `elj_cu()`: CUDA implementation;
- `elj_cl()`: OpenCL implementation, *if we would eventually add any OpenCL code
  in the parallel universe.*

Don’t put the suffix at the wrong place. Examples such as `empole_cu_nonewald()`
and `epolar_acc2_ewald()` should be corrected to `empole_nonewald_cu()` and
`epolar_ewald2_acc()`.

Variations of suffixes (e.g., `_acc1`, `_cu2`) are reserved for **kernels,
kernel templates, and function templates.** For examples, `edisp_acc1<Version,
PMEType>()` is the OpenACC kernel template used inside the regular function
`edisp_nonewald_acc()`; `edisp_cu1()` is the CUDA kernel template for the
dispersion energy. If you want to implement `func3()` and `func4()` in CUDA,
don’t name them `func_cu3()` and `func_cu4()`. Name them `func3_cu()` and
`func4_cu()` instead.
