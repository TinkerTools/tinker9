cmake_minimum_required (VERSION 3.12)


project (Tinker9_LibTinker LANGUAGES Fortran)


#'PATTERN' ==> [=[PATTERN]=]
set (__T9_GREP_ARGS [=[(?<=OBJS =)(.*\n)*.*(?=EXEFILES =)]=])
execute_process (
   COMMAND
      grep -Poz ${__T9_GREP_ARGS} Makefile
   WORKING_DIRECTORY
      ${PROJECT_SOURCE_DIR}/../tinker/make
   OUTPUT_VARIABLE
      __T9_LIBTINKER_O
)


string (REPLACE "\t" "" __T9_LIBTINKER_O ${__T9_LIBTINKER_O})
string (REPLACE "\n" "" __T9_LIBTINKER_O ${__T9_LIBTINKER_O})
string (REPLACE "\r" "" __T9_LIBTINKER_O ${__T9_LIBTINKER_O})
string (REPLACE ".o" ".f" __T9_LIBTINKER_O ${__T9_LIBTINKER_O})
string (REPLACE " \\" ";" __T9_LIBTINKER_O ${__T9_LIBTINKER_O})


# libtinker.a
foreach (var ${__T9_LIBTINKER_O})
   list (APPEND __T9_LIBTINKER_F ${PROJECT_SOURCE_DIR}/../tinker/source/${var})
endforeach ()
add_library (tinker STATIC ${__T9_LIBTINKER_F})
install (TARGETS tinker DESTINATION "${CMAKE_INSTALL_PREFIX}")
