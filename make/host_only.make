ifdef host
host := $(host)
else
host := 1
endif

ifeq ($(host),$(filter $(host),false off 0))
	host_only__ := false
	host_only_macro__ :=
else ifeq ($(host),$(filter $(host),true on 1))
	host_only__ := true
	host_only_macro__ := -DTINKER_HOST
endif
