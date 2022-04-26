#pragma once

//note that this file does not dispatch code to device or host version,
//but only facilitates the itellicense of device code

#if defined(__CUDACC__) || defined(__HIPCC__)
#define LEBEDEV_GROUP_INFO_BOTH_CALLABLE __host__ __device__
#define LEBEDEV_GROUP_INFO_GLOBAL __global__
#define LEBEDEV_GROUP_INFO_BOTH_INLINE __forceinline__
#else
#define LEBEDEV_GROUP_INFO_BOTH_CALLABLE
#define LEBEDEV_GROUP_INFO_GLOBAL
#define LEBEDEV_GROUP_INFO_BOTH_INLINE inline
#endif