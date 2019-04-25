############################################################
#                                                          #
#  Copy this Makefile and inc file to the build directory  #
#                                                          #
############################################################

config := default__
ifeq ($(config),default__)
	config := gnu.inc
endif
include $(config)

f77_default__ := gfortran

############################################################

ifeq ($(F77),f77)
	F77 := $(f77_default__)
endif
found__ := false

use_gfortran__ := false
f77__ := $(shell echo $(F77) | cut -c 1-8)
ifeq ($(f77__),gfortran)
	use_gfortran__ := true
	found__ := true
	fortran_name_macro := TINKER_GFORTRAN
endif

############################################################

os__ := $(shell uname -s)
