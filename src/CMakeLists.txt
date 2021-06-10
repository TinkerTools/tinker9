cmake_minimum_required (VERSION 3.12)


set (__T9_CPP_SRC
analyze_x.cpp
bar_x.cpp
box.cpp
cflux.cpp
compilers.cpp
couple.cpp
dynamic_x.cpp
eangle.cpp
eangtor.cpp
ebond.cpp
echarge.cpp
echglj.cpp
echgtrn.cpp
edisp.cpp
egeom.cpp
eimprop.cpp
eimptor.cpp
elec.cpp
emplar.cpp
empole_chgpen.cpp
empole.cpp
energy_buffer.cpp
energy.cpp
eopbend.cpp
epitors.cpp
epolar_chgpen.cpp
epolar.cpp
erepel.cpp
error_nix.cpp
estrbnd.cpp
estrtor.cpp
etors.cpp
etortor.cpp
eurey.cpp
evalence.cpp
evdw.cpp
exec.cpp
fc.cpp
field_chgpen.cpp
field.cpp
group.cpp
induce.cpp
induce_donly.cpp
initial.cpp
io_fort_str.cpp
io_text.cpp
lf_lpiston.cpp
lpiston.cpp
mdegv.cpp
mdintg.cpp
mdintg_respa.cpp
mdintg_verlet.cpp
mdpq.cpp
mdpt.cpp
mdsave.cpp
mechanic2.cpp
mechanic.cpp
minimize_x.cpp
mod.cpp
molecule.cpp
moments.cpp
nblist.cpp
nose.cpp
osrw.cpp
platform.cpp
pme.cpp
pmestuf.cpp
potent.cpp
random.cpp
rattle.cpp
rc_man.cpp
spatial2.cpp
spatial.cpp
switch.cpp
testgrad_x.cpp
tinker_mod.cpp
tinker_rt.cpp
)
if (HOST)
   set (__T9_CPP_HOST_SRC
host/darray.cpp
host/execq.cpp
host/fft.cpp
host/gpu_card.cpp
   )
   list (APPEND __T9_CPP_SRC "${__T9_CPP_HOST_SRC}")
else ()
   set (__T9_CPP_CUDART_SRC
cudart/darray.cpp
cudart/error.cpp
cudart/execq.cpp
cudart/fft.cpp
cudart/gpu_card.cpp
cudart/thrust_cache.cpp
   )
   list (APPEND __T9_CPP_SRC "${__T9_CPP_CUDART_SRC}")
endif ()
add_library (__t9_cpp_o OBJECT "${__T9_CPP_SRC}")
set_target_properties (__t9_cpp_o PROPERTIES
   CXX_STANDARD
      ${T9_CPPSTD}
)
target_compile_definitions (__t9_cpp_o PRIVATE "${T9_DEFS}")
target_include_directories (__t9_cpp_o SYSTEM PRIVATE "${T9_SYS_INCPATH}")
target_include_directories (__t9_cpp_o PRIVATE "${T9_INCPATH}")


add_library (__t9_version_related_o OBJECT
   info_x.cc
   promo.cc
)
set_target_properties (__t9_version_related_o PROPERTIES
   CXX_STANDARD
      ${T9_CPPSTD}
)
if (__T9_GIT_SHA1)
   target_compile_definitions (__t9_version_related_o PRIVATE
      TINKER9_GIT_SHORT_HASH=${__T9_GIT_SHA1})
endif()
if (__T9_GIT_DATE)
   target_compile_definitions (__t9_version_related_o PRIVATE
      TINKER9_GIT_DATE="${__T9_GIT_DATE}")
endif ()
target_compile_definitions (__t9_version_related_o PRIVATE "${T9_DEFS}")
target_include_directories (__t9_version_related_o SYSTEM PRIVATE "${T9_SYS_INCPATH}")
target_include_directories (__t9_version_related_o PRIVATE "${T9_INCPATH}")


add_library (tinker9_cpp STATIC
   $<TARGET_OBJECTS:__t9_cpp_o>
   $<TARGET_OBJECTS:__t9_version_related_o>
)
