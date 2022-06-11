# printf "set (__T9_OBJC_SRC\n%s\n)\n" "$(LC_COLLATE=C ls -1a *.m)" > cmakesrc.txt
cmake_minimum_required (VERSION 3.15)


enable_language (C)


include (cmakesrc.txt)
add_library (tinker9_objc OBJECT "${__T9_OBJC_SRC}")
