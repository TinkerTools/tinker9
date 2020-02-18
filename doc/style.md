# Style Guide

Tinker is now offloading part of future developments to C++.
Some of the signature features of C++ include:
   - most of the C features and procedural programming,
   - object-oriented programming,
   - generic programming (template),
   - the standard template library (STL),

which in my opinion, are only coarsely related,
and are in fact 4 different yet similar domain-specific languages.
I've heard about a comment on Fortran,
"New Fortran is great but no one uses it."
Similarly, C++ is great but no one can master all of it, and more
importantly, it doesn't mean we should abuse its features in a project.

Even the use of some of the popular and well-supported "good" features
should be limited. Everyone (and I) loves `std::unique_ptr` because
it is *almost* equivalent to raw pointer with (usually) unchangeable ownership,
and it will automatically deallocate the memory at the end of its lifecycle.
It is therefore very appealing to use it to manage the global arrays in this
project. But how many of you (and me) know how to let it:
   - cooperate with `cudaFree(void*)` function to auto-deallocate GPU memory
   (which is less commonly seen but still straightforward),
   - cooperate with `cudaMalloc(void**, size_t)` nice and cleanly,
   - be passed to `OpenACC` and `CUDA` kernals elegantly as raw pointer?

Because of the flexibility and complexity of the language, restrictions on the
coding style are necessary. This style guide will follow the structure of the
[Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html),
and many words were directly adopted from the Google C++ Style Guide,
but it never means all of these rules apply to this Tinker project.


## Table of Contents

| Sections |              |
|----------|:-------------|
| [C++ Version](#cppvers) | |
| [Header Files](#head)   | [Self-contained Headers](#hd.selfcontained) &emsp; [Header Guard](#hd.guard) &emsp; [Inlined Functions](#hd.inlfunc) &emsp; [Names and Order of Includes](#hd.inlord) |
| [Naming](#naming)       | [Constant Names](#nm.const) &emsp; [Namespace Names](#nm.namespace) &emsp; [Enumerator Names](#nm.enum) &emsp; [Macro Names](#nm.macro) |
| [Comments](#comment)    | |
| [Formatting](#format)   | |


<a name='cppvers'></a>
## C++ Version
Currently, code should target C++11.

C++14 and C++17 are good, but the support of the toolchain to the newer
standards may be limited, and more likely, be limited by the combination
several problems, e.g., the Nvidia driver is old and cannot be updated because
the GPU is old or it is difficult to update cluster's OS, so CUDA is limited
to 10.0, which only supports C++14 but not C++17.


<!--  -->


<a name='header'></a>
## Header Files

<a name='hd.selfcontained'></a>
### Self-contained Headers
Header files should be self-contained (compile on their own) and end in `.h`.
Non-header files that are meant for inclusion should end in `.hh` and be used
sparingly.

**Bad Example**
```cpp
// a.h: self-contained.
struct A {};

// b.h: not self-contained; a.h should have been included.
int b(const A&);

// c.cpp: compiles; won't compile if a.h is not included.
#include "a.h"
#include "b.h"
int use_b(const A& a) {
   return b(a);
}
```

**Good Example**
```cpp
// a.h: self-contained.
struct A {};

// b.h: self-contained.
#include "a.h"
int b(const A&);

// c.cpp: good.
#include "b.h"
int use_b(const A& a) {
   return b(a);
}
```

<a name='hd.guard'></a>
### Header Guard
All header files should have `#pragma once` guards to prevent multiple
inclusion. Your compiler and/or machine must be **very** old if they don't
support this feature, and Tinker won't be able to run on your machine anyway.

Don't do this, unless you have to:
```cpp
#ifndef TINKER_UNIQUE_HEADER_GUARD_WITH_FILENAME_LIKE_THIS_H
#define TINKER_UNIQUE_HEADER_GUARD_WITH_FILENAME_LIKE_THIS_H
// Do we include the DIRECTORY_ name in a FILENAME_HEADER_GUARD?
// Do we add underbar after the HEADER_GUARD, like HEADER_GUARD_?
// How many _s? ZERO? ONE_? TWO__? THREE___?
// ...
// And they will become black holes if you update the filenames.
#endif
```

<a name='hd.inlfunc'></a>
### Inlined Functions
Yes, `inline` is a keyword in C/C++, but what it mainly does is to please the
linker, not even the compiler. If you really want to put a function
implementation into a header file, fine, use this keyword, or the linker will
be confused by several copies of the machine code of this function in different
object files. On the other hand, inlining a function doesn't always improve
the CPU performance, so there is no guarantee that the compiler will do so.
If you find inlining a CPU function is critical to the performance of the code,
there are other ways to force the compilers to do that.

**Definition**
You can declare functions in a way that allows the compiler to expand them
inline rather than calling them through the usual function call mechanism.

**Pros**
Inlining a function can generate more efficient object code, as long as the
inlined function is small.

**Cons**
Overuse of inlining can actually make programs slower. Depending on a function's
size, inlining it can cause the code size to increase (for long functions) or
decrease (for very short functions). Time to compile the code can be very
different.

As for CUDA, things are a bit different, where the compiler has promised to
aggressively inline every `__device__` kernel. One evidence is you don't get
to see these `__device__` kernels when you profile the CUDA code. But you don't
always have to use `inline` for the `__device__` kernel either. Two examples
are as follows.

**Need to Inline**
```cpp
// cu_a.h: used by other CUDA kernels in different files.
__device__ inline int a() { return 1; }

// b.cu
#include "cu_a.h"
__global__ void b(int* t) { *t += a(); }

// c.cu
#include "cu_a.h"
__global__ void c(int* t) { *t += a(); }
```

**No Need to Inline**
```cpp
// d.cu: or put all of the kernels in the same file.
__device__ int a() { return 1; }
__global__ void b(int* t) { *t += a(); }
__global__ void c(int* t) { *t += a(); }
```

<a name='hd.inlord'></a>
### Names and Order of Includes


<!--  -->


<a name='naming'></a>
## Naming

Rule of thumb:
  - readability first,
  - names should not be too long, so abbreviations are acceptable, as long as
    they are:
      - widely used (`PBC`),
      - or straightforward with no ambiguity (`config` for configuration),
      - otherwise, document it first before use (`GK` for AMOEBA
        Generalized Kirkwood Model)
  - type names should use either `snake_case` or `UpperCamelCase`,
      - light-weighted types prefer `snake_case`, e.g. `real` is a redefined
        type of `double` or `float`,
      - heavy-weighted types (e.g. complicated classes) should always use
        `UpperCamelCase`, 
  - constant flags should look like `CONSTANT_FLAGS`,
  - no requirements for local variables except for readability.

<a name='nm.const'></a>
### Constant Names

`constexpr` and `const` can decorate more than just variables. I will not
explain the differences if they are used with functions, parameters, etc. Just
never take them for granted. `constexpr` *usually* means a variable can be
determined/evaluated at compile-time, whereas `const` variable only guarantees
a variable won't change its value after initialization, although compilers
usually are clever enough to know whether or not a `const` variable is also
`constexpr`. But still, use `constexpr` whenever it is possible. 

Constant variables can either be declared globally or within the scope of a
function and they either look like constant flags or local variables, choose
their forms based on the use.

**Examples**
```cpp
constexpr int CU_PLTFM = 0x002; // used as a constant flag.

double kcal_to_kJ(double kcal) {
   constexpr double kJ_per_kcal = 4.184; // used as a local variable.
   return kcal * kJ_per_kcal;
}
```

<a name='nm.namespace'></a>
### Namespace Names

<a name='nm.enum'></a>
### Enumerator Names

They are `CONSTANT_FLAGS` to improve the readability of the code.

Their types and underlying types are usually unspecified. When the type needs
specified, use `UpperCamelCase`.

<a name='nm.macro'></a>
### Macro Names

Don't use macros unless you have to.

Most macros should look like a `CONSTANT_FLAG`, even if they work like
functions. Only few macros can look like `regular_functions` if everyone loves
how they look.

**Examples**
```cpp
// just an example; it should have been used as a constexpr variable
#define ROUNDED_PI 3.0
// single precision square root
#define REAL_SQRT(x) sqrtf(x)
```

<!--  -->


<a name='comments'></a>
## Comments


<!--  -->


<a name='format'></a>
## Formatting
