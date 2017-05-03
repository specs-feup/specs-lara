//
// Created by jsevcik on 3/2/16.
//

#ifndef ROUTING_EVEHICLESACCESS_H
#define ROUTING_EVEHICLESACCESS_H

#include <type_traits>

#include "../Utility/Definitions.h"

using Routing::byte;

namespace Routing {

    enum class EVehiclesAccess : byte {
        NoneVeh = 0,
        Regular = 1,
        PublicService = 2,
        Trucks = 4, // nakladni nad 3,5 t
        LightCommercialVehicles = 8, // nakladni do 3,5 t
        LongerHeavierVehicle = 16, // nakladni do 44 tun
        All = Regular | PublicService | Trucks | LightCommercialVehicles | LongerHeavierVehicle
    };
}

inline Routing::EVehiclesAccess operator|(Routing::EVehiclesAccess lhs, Routing::EVehiclesAccess rhs) {
    return (Routing::EVehiclesAccess) (static_cast<std::underlying_type<Routing::EVehiclesAccess>::type>(lhs) |
                                       static_cast<std::underlying_type<Routing::EVehiclesAccess>::type>(rhs));
}

inline Routing::EVehiclesAccess &operator|=(Routing::EVehiclesAccess &lhs, Routing::EVehiclesAccess rhs) {
    lhs = (Routing::EVehiclesAccess) (static_cast<std::underlying_type<Routing::EVehiclesAccess>::type>(lhs) |
                                      static_cast<std::underlying_type<Routing::EVehiclesAccess>::type>(rhs));
    return lhs;
}

inline Routing::EVehiclesAccess operator&(Routing::EVehiclesAccess lhs, Routing::EVehiclesAccess rhs) {
    return (Routing::EVehiclesAccess) (static_cast<std::underlying_type<Routing::EVehiclesAccess>::type>(lhs) &
                                       static_cast<std::underlying_type<Routing::EVehiclesAccess>::type>(rhs)
    );
}

#endif //ROUTING_EVEHICLESACCESS_H
