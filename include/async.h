#pragma once
#include "rc_man.h"
#include <cstring>


TINKER_NAMESPACE_BEGIN
class StreamSt;
typedef StreamSt* Stream;


void async_acc_queue_data(rc_op);
extern void* async_acc;


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
