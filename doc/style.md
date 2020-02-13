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
| [Naming](#naming)       | |
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
support this feature, and Tinker won't be able to run on that machine anyway.

Don't do this:
```cpp
#ifndef TINKER_UNIQUE_HEADER_GUARD_WITH_FILENAME_LIKE_THIS_H
#define TINKER_UNIQUE_HEADER_GUARD_WITH_FILENAME_LIKE_THIS_H
// Do we include the DIRECTORY_ name in a FILENAME_HEADER_GUARD?
// Do we add underbar after the HEADER_GUARD, like HEADER_GUARD_?
// How many _s? ZERO? ONE_? TWO__? THREE___?
// ...
// And they will become black holes if you update the file names.
#endif
```


<a name='hd.inlfunc'></a>
### Inlined Functions

<a name='hd.inlord'></a>
### Names and Order of Includes


<!--  -->


<a name='naming'></a>
## Naming


<!--  -->


<a name='comments'></a>
## Comments


<!--  -->


<a name='format'></a>
## Formatting
