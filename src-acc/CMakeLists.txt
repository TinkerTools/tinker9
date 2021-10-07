cmake_minimum_required (VERSION 3.12)


## libtinker9_acc.a
set (__T9_ACC_SRC
alterchg_acc.cpp
box_acc.cpp
compilers_acc.cpp
cudalib_acc.cpp
dcflux_acc.cpp
device_zero_acc.cpp
eangle_acc.cpp
eangtor_acc.cpp
ebond_acc.cpp
echarge_acc.cpp
echgtrn_acc.cpp
edisp_acc.cpp
egeom_acc.cpp
ehal_acc.cpp
eimprop_acc.cpp
eimptor_acc.cpp
elj_acc.cpp
empole_chgpen_acc.cpp
empole_ewald_acc.cpp
empole_nonewald_acc.cpp
empole_opbend_acc.cpp
eopbend_acc.cpp
epitors_acc.cpp
epolar_chgpen_acc.cpp
epolar_dot_acc.cpp
epolar_ewald_acc.cpp
epolar_nonewald_acc.cpp
erepel_acc.cpp
estrbnd_acc.cpp
estrtor_acc.cpp
etors_acc.cpp
etortor_acc.cpp
eurey_acc.cpp
field_chgpen_acc.cpp
field_ewald_acc.cpp
field_nonewald_acc.cpp
induce_acc.cpp
induce_donly_acc.cpp
lf_lpiston_acc.cpp
lpiston_acc.cpp
mathfunc_parallel_acc.cpp
mddebug_acc.cpp
mdegv_acc.cpp
mdintg_acc.cpp
mdpq_acc.cpp
mdpt_acc.cpp
minimize_acc.cpp
nblist_acc.cpp
osrw_acc.cpp
pmestuf_acc.cpp
rattle_acc.cpp
rotpole_acc.cpp
torque_acc.cpp
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
