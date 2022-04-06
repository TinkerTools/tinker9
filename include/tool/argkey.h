#pragma once
#include <cstring>
#include <string>
#include <vector>

namespace tinker {
void nextarg(size_t len, char* str, int& exist);

template <size_t Len>
void nextarg(char (&str)[Len], int& exist)
{
   nextarg(Len, str, exist);
}

template <class T1, class T2>
void getKV(std::string k, T1& v, T2 vIfKNotFound);

template <class T>
void getKV(std::string k, std::vector<T>& v);
}
