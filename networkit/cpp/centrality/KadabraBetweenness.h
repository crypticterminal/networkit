/*
 * KadabraBetweenness.h
 *
 * Created on: 18.07.2018
 * 		 Author: Eugenio Angriman, Alexander van der Grinten
 */

#ifndef KADABRA_H_
#define KADABRA_H_

#include "../auxiliary/PrioQueue.h"
#include "../base/Algorithm.h"
#include "../graph/Graph.h"

namespace NetworKit {

class Status {
public:
	Status(const count k);
	const count k;
	std::vector<node> top;
	std::vector<double> approxTop;
	std::vector<bool> finished;
	std::vector<double> bet;
	std::vector<double> errL;
	std::vector<double> errU;
	count nPairs;
};

class SpSampler {
public:
	SpSampler(const Graph &G, const count seed);
	std::vector<node> randomPath();

private:
	const Graph &G;
	const count n;
	Graph pred;
	std::vector<count> ballInd;
	std::vector<count> dist;
	std::vector<count> q;
	std::vector<count> nPaths;

	inline node randomNode() const;
	void backtrackPath(const node u, const node v, const node start,
	                   std::vector<node> &path);
	void removeSomeEdges(const std::vector<node> &vertices, const count length);
	void removeAllEdges();
	count getDegree(const Graph &graph, node y, bool useDegreeIn);
};

class KadabraBetweenness : public Algorithm {
public:
	KadabraBetweenness(const Graph &G, const count k = 1,
	                   const double delta = 0.1, const double err = 0.01,
	                   count unionSample = 0, const count startFactor = 100);
	void run() override;
	std::vector<std::pair<node, double>> ranking() const;
	std::vector<node> topkNodesList() const;
	std::vector<double> topkScoresList() const;
	count getNumberOfIterations() const;

protected:
	const Graph &G;
	const count k;
	const double delta;
	const double err;
	const count n;
	const count startFactor;
	const bool absolute;
	count unionSample;
	count nPairs;
	double delta_l_min_guess;
	double delta_u_min_guess;
	double omega;

	std::vector<node> topkNodes;
	std::vector<double> topkScores;
	Aux::PrioQueue<int64_t, node> top;

	std::vector<double> approx;
	std::vector<double> delta_l_guess;
	std::vector<double> delta_u_guess;

	const double balancingFactor = 0.001;

	void init();
	void computeDeltaGuess();
	void computeBetErr(Status *status, std::vector<double> &bet,
	                   std::vector<double> &errL,
	                   std::vector<double> &errU) const;
	void oneRound(SpSampler &sampler);
	bool computeFinished(Status *status) const;
	void getStatus(Status *status) const;
	double computeF(const double btilde, const count iterNum,
	                const double deltaL) const;
	double computeG(const double btilde, const count iterNum,
	                const double deltaU) const;
	void fillPQ();
	void fillResult(Status *status);
};

inline std::vector<node> KadabraBetweenness::topkNodesList() const {
	return topkNodes;
}

inline std::vector<double> KadabraBetweenness::topkScoresList() const {
	return topkScores;
}

inline void KadabraBetweenness::fillResult(Status *status) {
	computeFinished(status);
	topkScores.reserve(n);
	topkNodes.reserve(n);
	for (count i = 0; i < k; ++i) {
		topkScores.push_back(status->approxTop[i] / (double)status->nPairs);
		topkNodes.push_back(status->top[i]);
	}

	double betk = status->approxTop[k - 1] / (double)status->nPairs;
	double lbetk =
	    betk - computeF(betk, status->nPairs, delta_l_guess[status->top[k - 1]]);
	count pos = k + 1;
	count i;
	for (i = k; i < status->k; ++i) {
		topkScores.push_back(status->approxTop[i] / (double)status->nPairs);
		topkNodes.push_back(status->top[i]);
	}
}

inline std::vector<std::pair<node, double>>
KadabraBetweenness::ranking() const {
	std::vector<std::pair<node, double>> result(topkNodes.size());
#pragma omp parallel for
	for (omp_index i = 0; i < result.size(); ++i) {
		result[i] = std::make_pair(topkNodes[i], topkScores[i]);
	}
	return result;
}

inline count KadabraBetweenness::getNumberOfIterations() const {
	return nPairs;
}

} // namespace NetworKit

#endif /* ifndef KADABRA_H_ */
