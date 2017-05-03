//
// Created by jsevcik on 3/3/16.
//

#ifndef ROUTING_RESTRICTEDTURN_H
#define ROUTING_RESTRICTEDTURN_H

#include <ostream>

namespace Routing {

    struct RestrictedTurn {

        int fromNodeId;
        int viaNodeId;
        int toNodeId;

        RestrictedTurn() : fromNodeId(-1), viaNodeId(-1), toNodeId(-1) {}

        RestrictedTurn(int fromNodeId, int viaNodeId, int toNodeId) :
                fromNodeId(fromNodeId),
                viaNodeId(viaNodeId),
                toNodeId(toNodeId) {};

        RestrictedTurn(const RestrictedTurn &other) : fromNodeId(other.fromNodeId), viaNodeId(other.viaNodeId),
                                                      toNodeId(other.toNodeId) {}

        RestrictedTurn(RestrictedTurn &&other) : fromNodeId(-1), viaNodeId(-1), toNodeId(-1) {
            fromNodeId = other.fromNodeId;
            viaNodeId = other.viaNodeId;
            toNodeId = other.toNodeId;

            other.fromNodeId = -1;
            other.viaNodeId = -1;
            other.toNodeId = -1;
        }

        std::ostream &operator<<(std::ostream &os) {
            os << "From: " << this->fromNodeId << ", Via: " << this->viaNodeId << ", To: " << this->toNodeId;

            return os;
        }

        friend std::ostream &operator<<(std::ostream &os, const RestrictedTurn &turn) {
            os << "From: " << turn.fromNodeId << ", Via: " << turn.viaNodeId << ", To: " << turn.toNodeId;

            return os;
        }

        bool operator==(const RestrictedTurn &other) {
            return fromNodeId == other.fromNodeId && viaNodeId == other.viaNodeId && toNodeId == other.toNodeId;
        }

        RestrictedTurn &operator=(const RestrictedTurn &other) {
            fromNodeId = other.fromNodeId;
            viaNodeId = other.viaNodeId;
            toNodeId = other.toNodeId;

            return *this;
        }
    };
}

#endif //ROUTING_RESTRICTEDTURN_H
