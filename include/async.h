#pragma once
#include "rc_man.h"
#include <cstring>


TINKER_NAMESPACE_BEGIN
/// \brief Deallocate the asynchronous stream.
void deallocate_stream(void*);
/// \brief Allocate the asynchronous stream.
void allocate_stream(void**);
/// \brief Synchronize the asynchronous stream.
void synchronize_stream(void*);


/// \brief
/// Copy between two device addresses without blocking the calling thread.
void copy_bytes_async(void* dst, const void* src, size_t nbytes, void* s);
TINKER_NAMESPACE_END
