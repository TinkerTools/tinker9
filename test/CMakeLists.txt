# LC_COLLATE=C ls -1a *.cpp
cmake_minimum_required (VERSION 3.12)


## test files
set (T9_TEST_SRC
angle.cpp
angtor.cpp
async.cpp
bond.cpp
bounds.cpp
box.cpp
chglj.cpp
chgtrn.cpp
disp.cpp
emhippo.cpp
ephippo.cpp
geom.cpp
improp.cpp
imptor.cpp
info.cpp
kinetic.cpp
local_frame.cpp
local_frame2.cpp
local_frame3.cpp
main.cpp
mathfunc.cpp
nacl.cpp
nblist.cpp
opbend.cpp
osrw.cpp
partial_charge.cpp
pitors.cpp
reduce.cpp
rephippo.cpp
strbnd.cpp
strtor.cpp
text.cpp
tinker_rt.cpp
torsion.cpp
tortor.cpp
trimatexp.cpp
truncated8.cpp
urey.cpp
vdw14.cpp
verlet.cpp
)
add_library (__t9_all_tests_o OBJECT ${T9_TEST_SRC})
target_include_directories (__t9_all_tests_o SYSTEM PRIVATE "${PROJECT_SOURCE_DIR}/ext/ext")
target_include_directories (__t9_all_tests_o PRIVATE        "${PROJECT_SOURCE_DIR}/include/syntax/acc")
target_link_libraries (__t9_all_tests_o PUBLIC tinker9_cpp)
