//
// Created by jsevcik on 3/2/16.
//

#ifndef ROUTING_ESPECIFICINFO_H
#define ROUTING_ESPECIFICINFO_H

#include <type_traits>

#include "../Utility/Definitions.h"

using Routing::ushort;

namespace Routing {
    enum ESpecificInfo : ushort {
        None = 0,
        TollWay = 1,
        SlipRoad = 2,
        MultiCarriageway = 4,
        Roundabout = 8,
        BorderCross = 16,
        Bridge = 32,
        Tunnel = 64
    };
}

inline Routing::ESpecificInfo operator|(Routing::ESpecificInfo lhs, Routing::ESpecificInfo rhs) {
    return (Routing::ESpecificInfo) (static_cast<std::underlying_type<Routing::ESpecificInfo>::type>(lhs) |
                                     static_cast<std::underlying_type<Routing::ESpecificInfo>::type>(rhs));
}

inline Routing::ESpecificInfo &operator|=(Routing::ESpecificInfo &lhs, Routing::ESpecificInfo rhs) {
    lhs = (Routing::ESpecificInfo) (static_cast<std::underlying_type<Routing::ESpecificInfo>::type>(lhs) |
                                    static_cast<std::underlying_type<Routing::ESpecificInfo>::type>(rhs));
    return lhs;
}

inline Routing::ESpecificInfo operator&(Routing::ESpecificInfo lhs, Routing::ESpecificInfo rhs) {
    return (Routing::ESpecificInfo) (static_cast<std::underlying_type<Routing::ESpecificInfo>::type>(lhs) &
                                     static_cast<std::underlying_type<Routing::ESpecificInfo>::type>(rhs));
}

#endif //ROUTING_ESPECIFICINFO_H
