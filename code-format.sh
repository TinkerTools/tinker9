#!/bin/bash

DETECTED=$(clang-format --version)            # Ubuntu clang-format version 10.0.1-++20200708123514+ef32c611aa2-1~exp1~20200707224111.189
DETECTED="${DETECTED##*clang-format version}" # 10.0.1-++20200708123514+ef32c611aa2-1~exp1~20200707224111.189
DETECTED="${DETECTED%%.*}"                    # 10
if [ $DETECTED -lt 10 ]; then
   echo Must use clang-format version 10.0.0+
   exit 1
fi

OS=$(uname -s)
if   [ $OS == Linux  ]; then
   SED='sed -i'
elif [ $OS == Darwin ]; then
   SED='perl -i -pe'
fi

PragmaDetected() {
   x=$(grep "#pragma" "$1")
   y=$(grep "#pragma once" "$1")
   if [ -z "$x" ]; then
      # no match
      return 1
   else
      nx=$(echo "$x" | wc -l)
      ny=$(echo "$y" | wc -l)
      if [ $ny -eq 1 ] && [ $nx -eq 1 ]; then
         # pragma once
         return 1
      else
         # found
         return 0
      fi
   fi
}

for x in "$@"; do

if PragmaDetected "$x"; then
   $SED 's/#pragma /\/\/#prag /g' "$x"
   clang-format -i -style=file "$x"
   $SED 's/\/\/ *#prag /#pragma /g' "$x"
else
   clang-format -i -style=file "$x"
fi

done

# clang-tidy
# ClangTidyFound=$(which clang-tidy)
# if [ "$ClangTidyFound" != "" ]; then
#    CMakeCacheFile=$(grep --include=\*.txt -irl "HOST:STRING=1" | head -1)
#    BuildDir="${CMakeCacheFile%%/*}"
#    for x in "$@"; do
#       FileExt="${x##*.}"
#       if [ "$FileExt" = "cpp" ] || [ "$FileExt" = "cxx" ] || [ "$FileExt" = "cc" ] ; then
#          clang-tidy -config= -header-filter=.* -p "$BuildDir" "$x"
#       fi
#    done
# fi
