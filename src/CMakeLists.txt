# LC_COLLATE=C ls -1a *.cpp {ff,ffpchg,math,md,tool}/*.cpp
cmake_minimum_required (VERSION 3.12)


set (__T9_CPP_SRC
cflux.cpp
compilers.cpp
echarge.cpp
echglj.cpp
echgtrn.cpp
edisp.cpp
elec.cpp
emplar.cpp
empole.cpp
empole_chgpen.cpp
epolar.cpp
epolar_chgpen.cpp
erepel.cpp
evdw.cpp
ff/atom.cpp
ff/box.cpp
ff/egvop.cpp
ff/elec.cpp
ff/energy.cpp
ff/molecule.cpp
ff/nblist.cpp
ff/pme.cpp
ff/potent.cpp
ff/spatial.cpp
ff/switch.cpp
ffpchg/eangle.cpp
ffpchg/eangtor.cpp
ffpchg/ebond.cpp
ffpchg/egeom.cpp
ffpchg/eimprop.cpp
ffpchg/eimptor.cpp
ffpchg/eopbend.cpp
ffpchg/epitors.cpp
ffpchg/estrbnd.cpp
ffpchg/estrtor.cpp
ffpchg/etors.cpp
ffpchg/etortor.cpp
ffpchg/eurey.cpp
ffpchg/evalence.cpp
field.cpp
field_chgpen.cpp
induce.cpp
induce_donly.cpp
initial.cpp
itgEnum.cpp
itgbBasic.cpp
itgbBerendsen.cpp
itgbLP22.cpp
itgbLogVAniso.cpp
itgbLogVIso.cpp
itgbMonteCarlo.cpp
itgbNhc06.cpp
itgiBasic.cpp
itgiLP22.cpp
itgiLeapFrogLP.cpp
itgiNhc06.cpp
itgiNhc96.cpp
itgiRespa.cpp
itgiVerlet.cpp
itgpBasic.cpp
itgpLogV.cpp
itgpRespa.cpp
itgtBasic.cpp
itgtBussi.cpp
itgtNhc.cpp
itgtNhc06.cpp
lf_lpiston.cpp
lpiston.cpp
math/lu.cpp
md/rattle.cpp
mdintg.cpp
mdppg.cpp
mdpq.cpp
mdpt.cpp
mdsave.cpp
mechanic2.cpp
mod.cpp
moments.cpp
nhc.cpp
nose.cpp
osrw.cpp
platform.cpp
random.cpp
rc_man.cpp
testrt.cpp
tinker_rt.cpp
tinker_suppl.cpp
tool/energybuffer.cpp
tool/errlinux.cpp
tool/errprint.cpp
tool/exec.cpp
tool/iofortstr.cpp
tool/iotext.cpp
trimatexp.cpp
xanalyze.cpp
xbar.cpp
xdynamic.cpp
xminimize.cpp
xtestgrad.cpp
)
if (HOST)
   list (APPEND __T9_CPP_SRC
host/darray.cpp
host/execq.cpp
host/fft.cpp
host/gpucard.cpp
   )
else ()
   list (APPEND __T9_CPP_SRC
cudart/darray.cpp
cudart/error.cpp
cudart/execq.cpp
cudart/fft.cpp
cudart/gpucard.cpp
cudart/thrustcache.cpp
   )
endif ()
add_library (tinker9_cpp OBJECT "${__T9_CPP_SRC}")


add_library (tinker9_version OBJECT
   xinfo.cc
   promo.cc
)
target_link_libraries (tinker9_version PUBLIC tinker9_cpp)


# macro
target_compile_definitions (tinker9_cpp PUBLIC ${T9_DEFS})
if (__T9_GIT_SHA1)
   target_compile_definitions (tinker9_version PRIVATE
      TINKER9_GIT_SHORT_HASH=${__T9_GIT_SHA1})
endif()
if (__T9_GIT_DATE)
   target_compile_definitions (tinker9_version PRIVATE
      TINKER9_GIT_DATE="${__T9_GIT_DATE}")
endif ()
# include
if (HOST)
else ()
   target_include_directories (tinker9_cpp PUBLIC "${CUDA_DIR}/include")
endif ()
target_include_directories (tinker9_cpp     PUBLIC  "${PROJECT_SOURCE_DIR}/include")
target_include_directories (tinker9_cpp     PRIVATE "${PROJECT_SOURCE_DIR}/include/syntax/acc")
target_include_directories (tinker9_version PRIVATE "${PROJECT_SOURCE_DIR}/include/syntax/acc")
# lib
target_link_libraries (tinker9_cpp PUBLIC tinkerFToCpp)
# FFTW
if (HOST)
   target_include_directories (tinker9_cpp SYSTEM PUBLIC "${FFTW_DIR}/include")
   target_link_libraries (tinker9_cpp PUBLIC t9_lfftw t9_lfftw_threads)
else ()
endif ()
