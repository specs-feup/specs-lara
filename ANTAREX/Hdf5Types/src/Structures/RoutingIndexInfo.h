#ifndef ROUTINGINDEXINFO_H_
#define ROUTINGINDEXINFO_H_

namespace Routing {

	struct IndexInfo {
		long long timeCreated;
		signed char isDirectMapped; //bool
		unsigned int numOfParts;
	};

	struct IndexPartInfo {
		int partNumber;
		float left;
		float top;
		float right;
		float bottom;
	};
}
#endif
