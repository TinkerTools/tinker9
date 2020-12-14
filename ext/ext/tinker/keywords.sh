#!/bin/bash

OUTPUT=keywords-fortran.txt
git -C ../../../tinker/ \
   log -1 --pretty=format:"# %h Tinker (Fortran)%n" \
   > "${OUTPUT}"

grep ' if (keyword(' ../../../tinker/source/*.f | while read line
do
   k="${line#* .eq. "'"}"
   k="${k#*.eq."'"}"
   k="${k%%" ')"*}"
   echo $k
done | sort | uniq >> "${OUTPUT}"

# while read line
# do
#    if [ ${line::1} != '#' ]; then
#       echo "**"$line"**"
#       echo
#       echo ".. index::" $line
#       echo
#    fi
# done < "${OUTPUT}"
