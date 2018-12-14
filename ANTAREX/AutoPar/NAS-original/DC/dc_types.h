#ifndef __DC_TYPES_H__
#define __DC_TYPES_H__


#ifdef WINNT
#ifndef HAS_INT64
typedef __int64             int64;
typedef int                 int32;
#endif
typedef unsigned __int64   uint64;
typedef unsigned int       uint32;
#else
#ifndef HAS_INT64
typedef long long           int64;
typedef int                 int32;
#endif
typedef unsigned long long uint64;
typedef unsigned int       uint32;
#endif

#endif