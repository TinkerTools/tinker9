cmake_minimum_required (VERSION 3.11)


file (GLOB __T9_CPP_SRC "*.cpp")
if (HOST)
   file (GLOB __T9_CPP_HOST_SRC "host/*.cpp")
   list (APPEND __T9_CPP_SRC "${__T9_CPP_HOST_SRC}")
else ()
   file (GLOB __T9_CPP_CUDART_SRC "cudart/*.cpp")
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
