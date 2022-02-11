cmake_minimum_required (VERSION 3.12)


## test files
set (T9_TEST_SRC
main_all_tests.cpp
test_angle.cpp
test_angtor.cpp
test_async.cpp
test_bond.cpp
test_bounds.cpp
test_box.cpp
test_chglj.cpp
test_chgtrn.cpp
test_disp.cpp
test_emhippo.cpp
test_ephippo.cpp
test_geom.cpp
test_improp.cpp
test_imptor.cpp
test_info.cpp
test_kinetic.cpp
test_local_frame2.cpp
test_local_frame3.cpp
test_local_frame.cpp
test_mathfunc.cpp
test_nacl.cpp
test_nblist.cpp
test_opbend.cpp
test_osrw.cpp
test_partial_charge.cpp
test_pitors.cpp
test_reduce.cpp
test_rephippo.cpp
test_rt.cpp
test_strbnd.cpp
test_strtor.cpp
test_text.cpp
test_tinker_rt.cpp
test_torsion.cpp
test_tortor.cpp
test_trimatexp.cpp
test_truncated8.cpp
test_urey.cpp
test_vdw14.cpp
test_verlet.cpp
)
add_library (__t9_all_tests_o OBJECT ${T9_TEST_SRC})
target_include_directories (__t9_all_tests_o SYSTEM PRIVATE "${PROJECT_SOURCE_DIR}/ext/ext")
target_include_directories (__t9_all_tests_o PRIVATE        "${PROJECT_SOURCE_DIR}/include/syntax/acc")
target_link_libraries (__t9_all_tests_o PUBLIC tinker9_cpp)
