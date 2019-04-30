fortran_name_macro := TINKER_GFORTRAN

# Fortran runtime
# ubuntu 16.04.6, apt install, gforran 5.4.0
fortrtflag := -L/usr/lib/gcc/x86_64-linux-gnu/5.4.0 -lgfortran -lgomp
