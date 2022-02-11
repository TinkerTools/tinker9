cmake_minimum_required (VERSION 3.12)


## libtinker9_acc.a
add_library (tinker9_acc OBJECT
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
target_include_directories (tinker9_acc PRIVATE "${PROJECT_SOURCE_DIR}/include/syntax/acc")
target_link_libraries (tinker9_acc PRIVATE tinker9_cpp)


########################################################################


if (HOST)
else ()
   string (APPEND CMAKE_CXX_FLAGS " CUDA_HOME=${CUDA_DIR}")
   target_compile_options (tinker9_acc PUBLIC
      -acc verystrict
      -Minfo=accel
   )
   ## Debug add flag: -ta=tesla:lineinfo,cc60,cc70
   ## Release add flag: -ta=tesla:fastmath,cc60,cc70
   target_compile_options (tinker9_acc PUBLIC
      "$<$<CONFIG:DEBUG>:-ta=tesla:lineinfo${__T9_ACC_CCLST4}>"
      "$<$<CONFIG:RELWITHDEBINFO>:-ta=tesla:lineinfo,fastmath${__T9_ACC_CCLST4}>"
      "$<$<CONFIG:RELEASE>:-ta=tesla:fastmath${__T9_ACC_CCLST4}>"
      "$<$<CONFIG:MINSIZEREL>:-ta=tesla:fastmath${__T9_ACC_CCLST4}>"
   )
endif ()
