ifdef fftw_lib
fftw_lib__ = $(fftw_lib)
else
ifdef fftw_dir
fftw_lib__ = $(fftw_dir)/lib
else
$(error fftw_dir is not set)
endif
endif
ifdef fftw_include
	fftw_include__ = $(fftw_include)
else
	ifdef fftw_dir
		fftw_include__ = $(fftw_dir)/include
	else
		$(error fftw_dir is not set)
	endif
endif
ifdef fftw_dir
	ifndef fftw_include__
		fftw_include__ = $(fftw_dir)/include
	endif
	ifndef fftw_lib__
		fftw_lib__ = $(fftw_dir)/lib
	endif
endif

