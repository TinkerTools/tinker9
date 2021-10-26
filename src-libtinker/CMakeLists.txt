cmake_minimum_required (VERSION 3.12)


#'PATTERN' ==> [=[PATTERN]=]
# The old command
# grep -Poz '(?<=OBJS =)(.*\n)*.*(?=EXEFILES =)' Makefile
# does not work with the BSD grep.
# sed -ne '/OBJS =/,/EXEFILES =/{s/OBJS =//; /EXEFILES =/!p;}' Makefile
set (__T9_GREP_ARGS [=[/OBJS =/,/EXEFILES =/{s/OBJS =//\; /EXEFILES =/!p\;}]=])
execute_process (
   COMMAND
      sed -ne ${__T9_GREP_ARGS} Makefile
   WORKING_DIRECTORY
      ${CMAKE_CURRENT_SOURCE_DIR}/../tinker/make
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
   list (APPEND __T9_LIBTINKER_F ${CMAKE_CURRENT_SOURCE_DIR}/../tinker/source/${var})
endforeach ()
add_library (tinker9_libtinker STATIC ${__T9_LIBTINKER_F})
