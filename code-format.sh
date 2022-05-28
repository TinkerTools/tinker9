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
   local x=$(grep "#pragma" "$1")
   if [ -z "$x" ]; then
      # if x is empty: no match
      return 1
   else
      # if x is not empty
      local nx=$(echo "$x" | wc -l)
      if [ $nx -eq 1 ] && grep -q "#pragma once" "$1"; then
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
   # copy x to y
   y="$x"-save-for-clang-format
   cp "$x" "$y"
   # format y
   $SED 's/#pragma /\/\/#prag /g' "$y"
   clang-format -i -style=file "$y"
   $SED 's/\/\/ *#prag /#pragma /g' "$y"
   # compare x and y
   if cmp --silent "$x" "$y"; then
      # x and y are identical
      rm "$y"
   else
      # x and y are different
      mv "$y" "$x"
   fi
else
   clang-format -i -style=file "$x"
fi
done
