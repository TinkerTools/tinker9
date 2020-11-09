#error This header file should never have been included.


/// \defgroup macro  Predefined Macros
/// \defgroup prec  Precisions of MM and MD Types
/// \defgroup acc_syntax  OpenACC Specific Code
/// \defgroup cuda_syntax  CUDA Specific Code
/// \defgroup io  I/O and Text
/// \defgroup error  Errors and Exceptions
/// \defgroup test  Unit Tests

/// \defgroup async  Asynchronous Flow Control
/// \defgroup rc  Resource: Pointer, Allocation, Deallocation, Queue

/// \defgroup math  Math
/// \defgroup parallel_algo  Parallel Algorithms
/// \defgroup spatial2  Spatial Decomposition Version 2


/// \defgroup mdegv  Energy, Gradient, and Virial Tensor


/*
/// \defgroup spatial  Spatial Decomposition
/// \defgroup geom  Geometrical Restraints
/// \defgroup vdw  Van der Waals (VDW)
/// \defgroup charge  Partial Charge Electrostatics
/// \defgroup mpole  Multipole Electrostatics
/// \defgroup polar  AMOEBA Polarization Electrostatics
/// \defgroup pme  Particle Mesh Ewald


/// \defgroup mdcalc  Control Bits for MM and MD Calculations
/// \defgroup mdpq  Atom Number, Momentum (p), and Coordinates (q)
/// \defgroup mdsave  Saving MD Trajectory Snapshots
/// \defgroup mdpt  Thermostats (T) and Barostats (P)
/// \defgroup mdintg  Integrators
/// \defgroup box  Periodic Boundary Box
/// \defgroup osrw  OSRW


/// \defgroup rand  Random Number
/// \defgroup nvidia  Nvidia GPU
/// \defgroup bindc    Fortran Tinker to C Interface
*/


//====================================================================//


/// \def TINKER_ICPC
/// \ingroup macro
/// Is defined when Intel C++ compiler was detected.
#define TINKER_ICPC
/// \def TINKER_GCC
/// \ingroup macro
/// Is defined when GNU C++ compiler was detected.
#define TINKER_GCC
/// \def TINKER_CLANG
/// \ingroup macro
/// Is defined when Clang C++ (either Apple or LLVM) compiler was detected.
#define TINKER_CLANG
/// \def TINKER_APPLE_CLANG
/// \ingroup macro
/// Is defined when Clang C++ (Apple) compiler was detected.
#define TINKER_APPLE_CLANG
/// \def TINKER_LLVM_CLANG
/// \ingroup macro
/// Is defined when Clang C++ (LLVM) compiler was detected.
#define TINKER_LLVM_CLANG
/* documented in macro.h */
#define TINKER_EXTERN_DEFINITION_FILE
