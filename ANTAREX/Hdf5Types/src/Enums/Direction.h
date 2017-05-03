//
// Created by jsevcik on 3/2/16.
//

#ifndef ROUTING_DIRECTION_H
#define ROUTING_DIRECTION_H

#include "../Utility/Definitions.h"

using Routing::byte;

namespace Routing {
    enum class Direction : byte {
        In = 0,
        Both = 1,
        Out = 2
    };
}

#endif //ROUTING_DIRECTION_H
