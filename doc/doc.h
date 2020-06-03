#error This header file should never have been included.


/// \defgroup style  Style Guide


/* sed '/PATTERN/ r doc/style.md' doc/doc.h > output */
/**
 * \ingroup style
 * \page style_md_pasted_here
 */


/// \defgroup macro  Predefined Macros
/// \defgroup acc_syntax  OpenACC Specific Code
/// \defgroup cuda_syntax  CUDA Specific Code


/// \defgroup geom  Geometrical Restraints
/// \defgroup vdw  Van der Waals (VDW)
/// \defgroup charge  Partial Charge Electrostatics
/// \defgroup mpole  Multipole Electrostatics
/// \defgroup polar  AMOEBA Polarization Electrostatics
/// \defgroup pme  Particle Mesh Ewald


/// \defgroup prec  Fine-Grained Definitions of MM and MD Types
/// \defgroup mdcalc  Control Bits for MM and MD Calculations
/// \defgroup mdpq  Atom Number, Momentum (p), and Coordinates (q)
/// \defgroup mdegv  Energy, Gradient, and Virial Tensor
/// \defgroup mdsave  Saving MD Trajectory Snapshots
/// \defgroup mdpt  Thermostats (T) and Barostats (P)
/// \defgroup mdintg  Integrators
/// \defgroup box  Periodic Boundary Box
/// \defgroup spatial  Spatial Decomposition
/// \defgroup osrw  OSRW


/// \defgroup math  Math
/// \defgroup rand  Random Number
/// \defgroup io  I/O and Text
/// \defgroup error  Errors and Exceptions
/// \defgroup mem  Memory and Pointers
/// \defgroup nvidia  Nvidia GPU
/// \defgroup bindc    Fortran Tinker to C Interface


/// \defgroup test  Unit Tests


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
