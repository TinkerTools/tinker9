ifeq ($(prec),$(filter $(prec),real4 4 32 float single s))
	prec__ := real4__
	precision_macro__ := -DTINKER_SINGLE_PRECISION
else ifeq ($(prec),$(filter $(prec),real8 8 64 double d))
	prec__ := real8__
	precision_macro__ := -DTINKER_DOUBLE_PRECISION
endif
