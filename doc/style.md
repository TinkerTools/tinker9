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
and a lot of words were directly copied from the Google C++ Style Guide,
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

C++14 and C++17 are good, But the support to the newer standards
may be restricted by the toolchain, or by the combination of the hardware
and toolchain, e.g. CUDA 10.0 only supports C++14 but not C++17.


<!--  -->


<a name='header'></a>
## Header Files

<a name='hd.selfcontained'></a>
### Self-contained Headers

<a name='hd.guard'></a>
### Header Guard

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
