//
// Created by jsevcik on 3/9/16.
//

#ifndef ROUTING_SHORTCUTPOSITION_H
#define ROUTING_SHORTCUTPOSITION_H

#include <ostream>

namespace Routing {

    struct ShortcutPosition {
    public:
        int edgeId;
        int position;
        int count;

        ShortcutPosition(int edgeId, int position, int count) :
                edgeId(edgeId),
                position(position),
                count(count) { };

        friend std::ostream &operator<<(std::ostream &os, Routing::ShortcutPosition &pos) {
            os << "ID: " << pos.edgeId << ", Position: " << pos.position << ", Count: " << pos.count;
            return os;
        }
    };
}

#endif //ROUTING_SHORTCUTPOSITION_H
