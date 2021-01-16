cmake_minimum_required (VERSION 3.11)


file (GLOB T9_CPP_SRC "*.cpp")
if (HOST)
   file (GLOB T9_CPP_HOST_SRC "host/*.cpp")
   list (APPEND T9_CPP_SRC "${T9_CPP_HOST_SRC}")
else ()
   file (GLOB TINKER9_CPP_CUDART_SRC "cudart/*.cpp")
   list (APPEND T9_CPP_SRC "${TINKER9_CPP_CUDART_SRC}")
endif ()
add_library (tinker9_cpp STATIC "${T9_CPP_SRC}")
set_target_properties (tinker9_cpp PROPERTIES
   CXX_STANDARD
      11
)
target_compile_definitions (tinker9_cpp PRIVATE "${T9_DEFS}")
target_include_directories (tinker9_cpp SYSTEM PRIVATE "${T9_SYS_INCPATH}")
target_include_directories (tinker9_cpp PRIVATE "${T9_INCPATH}")
