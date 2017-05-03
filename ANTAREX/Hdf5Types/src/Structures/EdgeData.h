//
// Created by jsevcik on 3/2/16.
//

#ifndef ROUTING_EDGEDATA_H
#define ROUTING_EDGEDATA_H

#include <limits>
#include <ostream>

#include "../Enums/EVehiclesAccess.h"
#include "../Utility/Definitions.h"
#include "../Enums/ESpecificInfo.h"

using Routing::byte;
using Routing::sbyte;
using Routing::ushort;

namespace Routing {

    class EdgeData {

    public:
        byte speed;
        byte funcClass;
        byte lanes;
        sbyte incline;
        EVehiclesAccess vehicleAccess;
        ESpecificInfo specificInfo;

        EdgeData(byte speed, byte funcClass, byte lanes, EVehiclesAccess vehAccess, ESpecificInfo specInfo) :
                speed(speed),
                funcClass(funcClass),
                lanes(lanes),
                incline(0),
                vehicleAccess(vehAccess),
                specificInfo(specInfo) {};

        EdgeData(byte speed, byte funcClass, byte lanes, EVehiclesAccess vehAccess, ESpecificInfo specInfo,
                 ushort maxWeight, ushort maxHeight, byte maxAxleLoad, byte maxWidth, byte maxLength, sbyte incline) :

                EdgeData(speed, funcClass, lanes, vehAccess, specInfo) {
            this->maxWeight = maxWeight;
            this->maxHeight = maxHeight;
            this->maxAxleLoad = maxAxleLoad;
            this->maxWidth = maxWidth;
            this->maxLength = maxLength;
            this->incline = incline;
        };

        EdgeData(byte speed, byte funcClass, byte lanes, EVehiclesAccess vehAccess, ESpecificInfo specInfo,
                 float maxWeight, float maxHeight, float maxAxleLoad, float maxWidth, float maxLength, sbyte incline) :
                EdgeData(speed, funcClass, lanes, vehAccess, specInfo) {
            this->SetMaxWeight(maxWeight);
            this->SetMaxHeight(maxHeight);
            this->SetMaxAxleLoad(maxAxleLoad);
            this->SetMaxWidth(maxWidth);
            this->SetMaxLenght(maxLength);
            this->incline = incline;
        };


        EdgeData(const EdgeData &other) :
                speed(other.speed),
                funcClass(other.funcClass),
                lanes(other.lanes),
                incline(other.incline),
                vehicleAccess(other.vehicleAccess),
                specificInfo(other.specificInfo),
                maxWeight(other.maxWeight),
                maxHeight(other.maxHeight),
                maxAxleLoad(other.maxAxleLoad),
                maxWidth(other.maxWidth),
                maxLength(other.maxLength) {};

        EdgeData(EdgeData &&other) {
            this->speed = other.speed;
            this->funcClass = other.funcClass;
            this->lanes = other.lanes;
            this->incline = other.incline;
            this->vehicleAccess = other.vehicleAccess;
            this->specificInfo = other.specificInfo;
            this->maxWeight = other.maxWeight;
            this->maxHeight = other.maxHeight;
            this->maxAxleLoad = other.maxAxleLoad;
            this->maxWidth = other.maxWidth;
            this->maxLength = other.maxLength;
        };

        friend std::ostream &operator<<(std::ostream &os, const EdgeData &edgeData) {
            os << "S: " << int(edgeData.speed) << ", F: " << int(edgeData.funcClass) << ", L: "
               << int(edgeData.lanes) << ", A: " << int(edgeData.vehicleAccess) << ", I: "
               << int(edgeData.specificInfo);

            return os;
        }

        inline float GetMaxWeight(void) const { return this->maxWeight / DIM_CONS_GET; };

        inline ushort GetMaxWeightRaw(void) const { return this->maxWeight; };

        void SetMaxWeight(float value) { this->maxWeight = static_cast<ushort>(value / DIM_CONS_SET); };

        inline float GetMaxHeight(void) const { return this->maxHeight / DIM_CONS_GET; };

        inline ushort GetMaxHeightRaw(void) const { return this->maxHeight; };

        void SetMaxHeight(float value) { this->maxHeight = static_cast<ushort>(value / DIM_CONS_SET); };

        inline float GetMaxAxleLoad(void) const { return this->maxAxleLoad / DIM_CONS_GET; };

        inline byte GetMaxAxleLoadRaw(void) const { return this->maxAxleLoad; };

        void SetMaxAxleLoad(float value) { this->maxAxleLoad = static_cast<byte>(value / DIM_CONS_SET); };

        inline float GetMaxWidth(void) const { return this->maxWidth / DIM_CONS_GET; };

        inline byte GetMaxWidthRaw(void) const { return this->maxWidth; };

        void SetMaxWidth(float value) { this->maxWidth = static_cast<byte>(value / DIM_CONS_SET); };

        inline float GetMaxLenght(void) const { return this->maxLength / DIM_CONS_GET; };

        inline byte GetMaxLenghtRaw(void) const { return this->maxLength; };

        void SetMaxLenght(float value) { this->maxLength = static_cast<byte>(value / DIM_CONS_SET); };

    private:
        ushort maxWeight = std::numeric_limits<ushort>::max();
        ushort maxHeight = std::numeric_limits<ushort>::max();
        byte maxAxleLoad = std::numeric_limits<byte>::max();
        byte maxWidth = std::numeric_limits<byte>::max();
        byte maxLength = std::numeric_limits<byte>::max();

        const float DIM_CONS_GET = 10.0;
        const int DIM_CONS_SET = 10;
    };
}

#endif //ROUTING_EDGEDATA_H
