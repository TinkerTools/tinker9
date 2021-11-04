#pragma once

#if !defined(TINKER_IFORT) && !defined(TINKER_GFORTRAN)
#   warning We assume you used gfortran to compile the Tinker library. \
If so, add the compiler flag -DTINKER_GFORTRAN to mute this message, \
or -DTINKER_IFORT if the Intel Fortran compiler is used. Otherwise, \
you should implement the name mangling macro (i.e., TINKER_MOD) here.
#   define TINKER_GFORTRAN
#endif

#if defined(TINKER_GFORTRAN)
#   define TINKER_MOD(mod, var) __##mod##_MOD_##var
#elif defined(TINKER_IFORT)
#   define TINKER_MOD(mod, var) mod##_mp_##var##_
#endif
