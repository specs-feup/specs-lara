//
// Created by jsevcik on 2/23/16.
//

#ifndef ROUTING_NODERAW_H
#define ROUTING_NODERAW_H

#include <ostream>

#include "../Utility/Definitions.h"

using Routing::byte;

namespace Routing {

    struct NodeRaw final {
    public:
        int id;
        int latitudeInt;
        int longitudeInt;
        byte edgesCount;
        int edgesIndex;

        NodeRaw(int ID, int LatitudeInt, int LongitudeInt) {
            id = ID;
            latitudeInt = LatitudeInt;
            longitudeInt = LongitudeInt;
            edgesCount = 0;
            edgesIndex = -1;
        };

        NodeRaw(int ID, int LatitudeInt, int LongitudeInt, byte EdgesCount, int EdgesIndex) {
            id = ID;
            latitudeInt = LatitudeInt;
            longitudeInt = LongitudeInt;
            edgesCount = EdgesCount;
            edgesIndex = EdgesIndex;
        };

        NodeRaw(int ID, float Latitude, float Longitude) {
            id = ID;
            latitudeInt = ConvertFloatToE6Int(Latitude);
            longitudeInt = ConvertFloatToE6Int(Longitude);
            edgesCount = 0;
            edgesIndex = -1;
        };


        NodeRaw(int ID, float Latitude, float Longitude, byte EdgesCount, int EdgesIndex) {
            id = ID;
            latitudeInt = ConvertFloatToE6Int(Latitude);
            longitudeInt = ConvertFloatToE6Int(Longitude);
            edgesCount = EdgesCount;
            edgesIndex = EdgesIndex;
        };

        NodeRaw(const NodeRaw &other) :
                id(other.id),
                latitudeInt(other.latitudeInt),
                longitudeInt(other.longitudeInt),
                edgesCount(other.edgesCount),
                edgesIndex(other.edgesIndex) { };

        NodeRaw &operator=(const NodeRaw &other) {
            this->id = other.id;
            this->latitudeInt = other.latitudeInt;
            this->longitudeInt = other.longitudeInt;
            this->edgesCount = other.edgesCount;
            this->edgesIndex = other.edgesIndex;

            return *this;
        }

        friend std::ostream &operator<<(std::ostream &os, const NodeRaw &node) {
            os << "ID: " << node.id << ", Lng: " << node.GetLongitude() << ", Lat: " << node.GetLatitude() <<
            ", ECount: " << static_cast<int>(node.edgesCount) << ", EIndex: " << node.edgesIndex;
            return os;
        }

        float GetLatitude(void) const { return ConvertE6IntToFloat(this->latitudeInt); };

        void SetLatitude(const float value) {
            this->latitudeInt = static_cast<int32_t>(value * this->coordinateDividor);
        };

        float GetLongitude(void) const { return ConvertE6IntToFloat(this->longitudeInt); };

        void SetLongitude(const float value) {
            this->longitudeInt = static_cast<int32_t>(value * this->coordinateDividor);
        };

    private:
        const float coordinateDividor = 1E6f; // 10^6

        inline int ConvertFloatToE6Int(const float value) const {
            return static_cast<int32_t>(value * this->coordinateDividor);
        }

        inline float ConvertE6IntToFloat(const int value) const {
            return value / this->coordinateDividor;
        }
    };
}

#endif //ROUTING_NODERAW_H
