//
// Created by jsevcik on 3/2/16.
//

#ifndef ROUTING_INDEXGEOMETRYSHORTCUT_H
#define ROUTING_INDEXGEOMETRYSHORTCUT_H

#include <memory>
#include <vector>

#include "Geo/GpsPoint.h"

namespace Routing {

    struct IndexGeometryShortcut {
    public:
        int nodeId1;
        int nodeId2;

        IndexGeometryShortcut(int id1, int id2, const std::vector<float> &gpsPointsLatitude, const std::vector<float> &gpsPointsLongitude) :
                nodeId1(id1),
                nodeId2(id2),
                gpsPointsLatitude(new std::vector<float>(gpsPointsLatitude)),
				gpsPointsLongitude(new std::vector<float>(gpsPointsLongitude)) {

        };

        IndexGeometryShortcut(int &&id1, int &&id2, std::vector<float> &&gpsPointsLatitude, std::vector<float> &&gpsPointsLongitude) :
                nodeId1(std::move(id1)),
                nodeId2(std::move(id2)),
                gpsPointsLatitude(new std::vector<float>(std::move(gpsPointsLatitude))),
				gpsPointsLongitude(new std::vector<float>(std::move(gpsPointsLongitude)))	{ };

        std::ostream &operator<<(std::ostream &os) {
            os << "ID1: " << this->nodeId1 << ", ID2: " <<
            this->nodeId2 << ", Points: " << this->gpsPointsLatitude->size();

            return os;
        }

        std::vector<float> &GetGpsPointsLatitude(void) const { return *this->gpsPointsLatitude; };
		std::vector<float> &GetGpsPointsLongitude(void) const { return *this->gpsPointsLongitude; };

    private:
	    //std::unique_ptr<std::vector<GpsPoint>> points;
		//Changed field 'points' to conform with example in 'RoutingIndexInfoCompType.cpp'
		std::unique_ptr<std::vector<float>> gpsPointsLatitude;
		std::unique_ptr<std::vector<float>> gpsPointsLongitude;
    };

}

#endif //ROUTING_INDEXGEOMETRYSHORTCUT_H
