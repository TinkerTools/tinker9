#!/bin/bash

OS=$(uname -s)
if   [ $OS == Linux  ]; then
   SED='sed -i'
elif [ $OS == Darwin ]; then
   SED='perl -i -pe'
fi

for x in "$@"; do
   $SED 's/#pragma /\/\/#prag /g' "$x"
   clang-format -i -style=file "$x"
   $SED 's/\/\/ *#prag /#pragma /g' "$x"
done
