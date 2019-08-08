os__ := $(shell uname -s 2>/dev/null | tr "[:upper:]" "[:lower:]")
ifeq ($(os__),linux)
	shared_lib_suffix__ := so
else ifeq ($(os__),darwin)
	shared_lib_suffix__ := dylib
endif
