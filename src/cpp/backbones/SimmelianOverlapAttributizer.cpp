/*
 * SimmelianOverlapAttributizer.cpp
 *
 *  Created on: 22.07.2014
 *      Author: Gerd Lindner
 */

#include "SimmelianOverlapAttributizer.h"
#include <limits>

namespace NetworKit {

SimmelianOverlapAttributizer::SimmelianOverlapAttributizer(count maxRank) :
		maxRank(maxRank) {}

EdgeAttribute SimmelianOverlapAttributizer::getAttribute(const Graph& graph, const EdgeAttribute& attribute) {
	std::vector<RankedNeighbors> neighbors = getRankedNeighborhood(graph, attribute);
	EdgeAttribute overlapAttribute(graph.upperEdgeIdBound(), 0.0);

	graph.forEdges([&](node u, node v, edgeid eid) {
		Redundancy redundancy = getOverlap(u, v, neighbors, maxRank);

		overlapAttribute[eid] = std::max((double) redundancy.overlap, overlapAttribute[eid]);
	});

	return overlapAttribute;
}

} /* namespace NetworKit */
