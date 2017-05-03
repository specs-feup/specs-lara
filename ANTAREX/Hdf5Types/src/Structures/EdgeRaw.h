//
// Created by jsevcik on 3/9/16.
//

#ifndef ROUTING_EDGERAW_H
#define ROUTING_EDGERAW_H

#include <sys/types.h>
#include <ostream>

#include "../Enums/Direction.h"

namespace Routing {

    struct EdgeRaw {
        int edgeId;
        int nodeId;
        int length;
        Direction direction;
        ushort dataIndex;

        EdgeRaw(int edgeId, int nodeId, int length, Direction direction, ushort dataIndex) :
                edgeId(edgeId),
                nodeId(nodeId),
                length(length),
                direction(direction),
                dataIndex(dataIndex) { };

        friend std::ostream &operator<<(std::ostream &os, const EdgeRaw &edgeRaw) {
            os << "E: " << edgeRaw.edgeId << ", N: " << edgeRaw.nodeId << ", L: " << edgeRaw.length <<
            ", D: " << (std::underlying_type<Routing::Direction >::type)edgeRaw.direction << ", DI: " << edgeRaw.dataIndex;

            return os;
        };
    };
}

#endif //ROUTING_EDGERAW_H
