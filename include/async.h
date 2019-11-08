#pragma once
#include "macro.h"
#include <cstring>


TINKER_NAMESPACE_BEGIN
class StreamSt;
typedef StreamSt* Stream;


/// \brief Deallocate the asynchronous stream.
void deallocate_stream(Stream);
/// \brief Allocate the asynchronous stream.
void allocate_stream(Stream*);
/// \brief Synchronize the asynchronous stream.
void synchronize_stream(Stream);


/// \brief
/// Copy between two device addresses without blocking the calling thread.
void copy_bytes_async(void* dst, const void* src, size_t nbytes, Stream s);
TINKER_NAMESPACE_END
