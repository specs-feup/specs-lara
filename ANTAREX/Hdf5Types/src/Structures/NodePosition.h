//
// Created by jsevcik on 3/9/16.
//

#ifndef ROUTING_NODEPOSITION_H
#define ROUTING_NODEPOSITION_H

#include <ostream>

namespace Routing {

    struct NodePosition {
        int nodeId;
        int partNumber;
        int nodeIndex;

        NodePosition(int nodeId, int partNumber, int nodeIndex) :
                nodeId(nodeId),
                partNumber(partNumber),
                nodeIndex(nodeIndex) { };

        std::ostream &operator<<(std::ostream &os) {
            os << "ID: " << this->nodeId << ", Part: " << this->partNumber
            << ", Index: <<" << this->nodeIndex;

            return os;
        }
    };
}

#endif //ROUTING_NODEPOSITION_H
