#ifndef ROUTING_DEFINITIONS_H
#define ROUTING_DEFINITIONS_H

#include <sys/types.h>
#include <cstdint>

#if defined(__GNUC__) || defined(__ICL) || defined(__clang__)
#define expect(x, y) (__builtin_expect((x),(y)))
#else
#define expect(x, y) (x)
#endif

#define likely(x) (expect(x, 1))
#define unlikely(x) (expect(x, 0))

typedef unsigned char test;

namespace Routing {

    typedef unsigned char byte;
    typedef signed char sbyte;
    typedef uint16_t ushort;
}

#endif // ROUTING_DEFINITIONS_H
