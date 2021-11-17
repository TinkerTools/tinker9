#pragma once
#include "tool/io_fort_str.h"
#include "tool/io_print.h"
#include "tool/io_read.h"
#include "tool/io_text.h"
#include <tinker/routines.h>


namespace tinker {
void mechanic2();


void nextarg(size_t len, char* str, int& exist);


template <size_t Len>
void nextarg(char (&str)[Len], int& exist)
{
   nextarg(Len, str, exist);
}


template <class T1, class T2>
void get_kv(std::string k, T1& v, T2 vdefault);


template <class T>
void get_kbool(std::string k, T& v, bool v_if_k_not_found);
}
