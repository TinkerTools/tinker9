cmake_minimum_required (VERSION 3.12)


## libtinker9_acc.a
set (__T9_ACC_SRC
alterchg.cpp
box.cpp
compilers.cpp
cudalib.cpp
dcflux.cpp
device_zero.cpp
eangle.cpp
eangtor.cpp
ebond.cpp
echarge.cpp
echgtrn.cpp
edisp.cpp
egeom.cpp
ehal.cpp
eimprop.cpp
eimptor.cpp
elj.cpp
empole_chgpen.cpp
empole_ewald.cpp
empole_nonewald.cpp
empole_opbend.cpp
eopbend.cpp
epitors.cpp
epolar_chgpen.cpp
epolar_dot.cpp
epolar_ewald.cpp
epolar_nonewald.cpp
erepel.cpp
estrbnd.cpp
estrtor.cpp
etors.cpp
etortor.cpp
eurey.cpp
field_chgpen.cpp
field_ewald.cpp
field_nonewald.cpp
induce.cpp
induce_donly.cpp
lf_lpiston.cpp
lpiston.cpp
mathfunc_parallel.cpp
mddebug.cpp
mdegv.cpp
mdintg.cpp
mdpq.cpp
mdpt.cpp
minimize.cpp
nblist.cpp
osrw.cpp
pmestuf.cpp
rattle.cpp
rotpole.cpp
torque.cpp
)
add_library (tinker9_acc STATIC "${__T9_ACC_SRC}")
set_target_properties (tinker9_acc PROPERTIES
   CXX_STANDARD
      ${T9_CPPSTD}
)
target_compile_definitions (tinker9_acc PRIVATE "${T9_DEFS}")
target_include_directories (tinker9_acc SYSTEM PRIVATE "${T9_SYS_INCPATH}")
target_include_directories (tinker9_acc PRIVATE "${T9_INCPATH}")


########################################################################
if (HOST)
   if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      string (APPEND CMAKE_CXX_FLAGS_DEBUG " -Wall -Wno-unknown-pragmas")
   endif ()
else ()
   ## Compute Capability 60,70 -> ,cc60,cc70
   set (__T9_ACC_CCLST4) # ""
   foreach (var ${T9_CUCCLIST})
      string (APPEND __T9_ACC_CCLST4 ,cc${var})
   endforeach () # ,cc60,cc70


   target_compile_options (tinker9_acc PRIVATE
      CUDA_HOME=${CUDA_DIR}
      -acc verystrict
      -Minfo=accel
   )
   ## Debug add flag: -ta=tesla:lineinfo,cc60,cc70
   ## Release add flag: -ta=tesla:fastmath,cc60,cc70
   string (APPEND CMAKE_CXX_FLAGS_DEBUG          " -ta=tesla:lineinfo${__T9_ACC_CCLST4}")
   string (APPEND CMAKE_CXX_FLAGS_RELWITHDEBINFO " -ta=tesla:lineinfo,fastmath${__T9_ACC_CCLST4}")
   string (APPEND CMAKE_CXX_FLAGS_RELEASE        " -ta=tesla:fastmath${__T9_ACC_CCLST4}")
   string (APPEND CMAKE_CXX_FLAGS_MINSIZEREL     " -ta=tesla:fastmath${__T9_ACC_CCLST4}")
endif ()
