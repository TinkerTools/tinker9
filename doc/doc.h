#error This header file should never have been included.


/// \defgroup lang  Compilers and Syntax
/// \defgroup compiler  Compilers
/// \ingroup lang
/// \defgroup cpp_syntax  C++ Syntax
/// \ingroup lang
/// \defgroup cuda_syntax  CUDA Specific Code
/// \ingroup lang
/// \defgroup acc_syntax  OpenACC Specific Code
/// \ingroup lang


/// \defgroup prec  Precisions


/// \defgroup io  I/O and Text
/// \defgroup error  Errors and Exceptions
/// \defgroup test  Unit Tests
/// \defgroup async  Asynchronous Flow Control
/// \defgroup rc  Resource: Pointer, Allocation, Deallocation, Queue
/// \defgroup math  Math
/// \defgroup parallel_algo  Parallel Algorithms
/// \defgroup spatial2  Spatial Decomposition Version 2

/// \defgroup box  Periodic Boundary Box
/// \defgroup mdegv  Energy, Gradient, and Virial Tensor

/// \defgroup fft  Fast Fourier Transform

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
/// \defgroup osrw  OSRW


/// \defgroup rand  Random Number
/// \defgroup nvidia  NVIDIA GPU
*/


//====================================================================//


/// \ingroup compiler
/// \brief Macro for the Intel C++ compiler.
#define TINKER_ICPC
/// \ingroup compiler
/// \brief Macro for the GNU C++ compiler.
#define TINKER_GCC
/// \ingroup compiler
/// \brief Macro for the Clang C++ (either Apple or LLVM) compiler.
#define TINKER_CLANG
/// \ingroup compiler
/// \brief Macro for the Clang C++ (Apple) compiler.
#define TINKER_APPLE_CLANG
/// \ingroup compiler
/// \brief Macro for the Clang C++ (LLVM) compiler.
#define TINKER_LLVM_CLANG
/// \ingroup compiler
/// \brief Macro for the PGI or NVHPC C++ compiler.
#define TINKER_PGI
#define TINKER_EXTERN_DEFINITION_FILE
