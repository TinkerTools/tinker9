fortran_search_dir0__ := $(shell ($(fortran_compiler) -print-search-dirs | grep 'libraries: =' | sed -e 's/libraries: =//g' -e 's/:/ /g'))
fortran_search_dir1__ := $(shell (find $(fortran_search_dir0__) -type f -name libgfortran.a 2>/dev/null | head -1))
fortran_search_dir2__ := $(shell dirname $(fortran_search_dir1__))
fortran_search_dir3__ := $(shell (cd "$(fortran_search_dir2__)" && pwd -P))


shared_flags__ += -DTINKER_GFORTRAN
link_flags__ += -L$(fortran_search_dir3__) -lgfortran
