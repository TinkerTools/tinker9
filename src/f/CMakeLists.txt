cmake_minimum_required (VERSION 3.15)


add_library (tinker9_f OBJECT tinkersupplement.F)
target_link_libraries (tinker9_f PRIVATE tinker9_cpp)
