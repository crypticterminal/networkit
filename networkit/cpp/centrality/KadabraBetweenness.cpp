/*
 * KadabraBetweenness.cpp
 *
 * Created on: 18.07.2018
 * 		 Author: Eugenio Angriman, Alexander van der Grinten
 */

#include <cmath>
#include <omp.h>

#include "../auxiliary/Random.h"
#include "../auxiliary/Timer.h"
#include "../distance/APSP.h"
#include "../distance/Diameter.h"
#include "KadabraBetweenness.h"

namespace NetworKit {

Status::Status(const count k)
		: k(k), top(std::vector<node>(k)), approxTop(std::vector<double>(k)),
			finished(std::vector<bool>(k)), bet(std::vector<double>(k)),
			errL(std::vector<double>(k)), errU(std::vector<double>(k)) {}

KadabraBetweenness::KadabraBetweenness(const Graph &G, const count k,
																			 const double delta, const double err,
																			 count unionSample,
																			 const count startFactor)
		: G(G), k(k), delta(delta), err(err), n(G.upperNodeIdBound()),
			startFactor(startFactor), unionSample(unionSample), absolute(k == 0),
			top(n, n) {
	if (k > n) {
		throw std::runtime_error(
				"k is higher than the number of nodes of the input graph!");
	}

	if (delta >= 1 || delta <= 0) {
		throw std::runtime_error(
				"Delta should be greater than 0 and smaller than 1.");
	}

	if (err >= 1 || err <= 0) {
		throw std::runtime_error(
				"The error should be greater than 0 and smaller than 1.");
	}
}

bool KadabraBetweenness::computeFinished(Status *status) const {
	auto bet = status->bet;
	auto errL = status->errL;
	auto errU = status->errU;
	bool allFinished = true;

	count i;
	for (i = 0; i < status->k - 1; ++i) {
		bet[i] = status->approxTop[i] / (double)status->nPairs;
		errL[i] = computeF(bet[i], status->nPairs, delta_l_guess[status->top[i]]);
		errU[i] = computeG(bet[i], status->nPairs, delta_u_guess[status->top[i]]);
	}

	bet[i] = status->approxTop[i] / (double)status->nPairs;
	errL[i] = computeF(bet[i], status->nPairs, this->delta_l_min_guess);
	errU[i] = computeG(bet[i], status->nPairs, this->delta_u_min_guess);

	if (absolute) {
		for (count i = 0; i < status->k; ++i) {
			status->finished[i] = (errL[i] < err && errU[i] < err);
			allFinished = allFinished && status->finished[i];
		}
	} else {
		for (count i = 0; i < status->k; ++i) {
			if (i == 0) {
				status->finished[i] = (bet[i] - errL[i] > bet[i + 1] + errU[i + 1]);
			} else if (i < k) {
				status->finished[i] = (bet[i - 1] - errL[i - 1] > bet[i] + errU[i]) &&
															(bet[i] - errL[i] > bet[i + 1] + errU[i + 1]);
			} else {
				status->finished[i] = bet[k - 1] - errU[k - 1] > bet[i] + errU[i];
			}
			status->finished[i] =
					status->finished[i] || (errL[i] < err && errU[i] < err);
			allFinished = allFinished && status->finished[i];
		}
	}

	return allFinished;
}

// Computes the function f that bounds the betweenness of a vertex from below.
// For more information, see Borassi, Natale (2016).
double KadabraBetweenness::computeF(const double btilde, const count iterNum,
																		const double deltaL) const {
	double tmp = (((double)omega) / iterNum - 1. / 3);
	double errChern = (std::log(1. / deltaL)) * 1. / iterNum *
										(-tmp + std::sqrt(tmp * tmp + 2 * btilde * omega /
																											(std::log(1. / deltaL))));
	return std::min(errChern, btilde);
}

// Computes the function g that bounds the betweenness of a vertex from above.
// For more information, see Borassi, Natale (2016).
double KadabraBetweenness::computeG(const double btilde, const count iterNum,
																		const double deltaU) const {
	double tmp = (((double)omega) / iterNum + 1. / 3);
	double errChern = (std::log(1. / deltaU)) * 1. / iterNum *
										(tmp + std::sqrt(tmp * tmp + 2 * btilde * omega /
																										 (std::log(1. / deltaU))));
	return std::min(errChern, 1 - btilde);
}

void KadabraBetweenness::oneRound(SpSampler &sampler) {
	auto path = sampler.randomPath();
#pragma omp critical
	{
		++nPairs;
		for (node u : path) {
			++approx[u];
			top.insert(u, approx[u]);
		}
	}
}

void KadabraBetweenness::getStatus(Status *status) const {
#pragma omp critical
	{
		if (status != NULL) {
			for (count i = 0; i < unionSample; ++i) {
				status->top[i] = top.get_element(i);
				status->approxTop[i] = approx[status->top[i]];
			}
			status->nPairs = nPairs;
		}
	}
}

void KadabraBetweenness::computeBetErr(Status *status, std::vector<double> &bet,
																			 std::vector<double> &errL,
																			 std::vector<double> &errU) const {
	count i;
	double maxErr = std::sqrt(startFactor) * err / 4.;

	for (i = 0; i < k; ++i) {
		bet[i] = status->approxTop[i] / (double)status->nPairs;
	}

	if (absolute) {
		for (i = 0; i < status->k; ++i) {
			errL[i] = err;
			errU[i] = err;
		}
	} else {
		errU[0] = std::max(err, (bet[0] - bet[1]) / 2.);
		errL[0] = 10;
		for (i = 1; i < k; ++i) {
			errL[i] = std::max(err, (bet[i - 1] - bet[i]) / 2.);
			errU[i] = std::max(err, (bet[i] - bet[i + 1]) / 2.);
		}
		for (i = k; i < status->k; ++i) {
			errL[i] = 10;
			errU[i] = std::max(err, bet[k - 1] + (bet[k - 1] - bet[k]) / 2. - bet[i]);
		}
		for (i = 0; i < k - 1; ++i) {
			if (bet[i] - bet[i + 1] < maxErr) {
				errL[i] = err;
				errU[i] = err;
				errL[i + 1] = err;
				errU[i + 1] = err;
			}
		}
		for (i = k + 1; i < status->k; ++i) {
			if (bet[k] - bet[i] < maxErr) {
				errL[k] = err;
				errU[k] = err;
				errL[i] = err;
				errU[i] = err;
			}
		}
	}
}

void KadabraBetweenness::computeDeltaGuess() {
	double a = 0,
				 b = 1. / err / err * std::log(n * 4 * (1 - balancingFactor) / delta),
				 c = (a + b) / 2;
	double sum;

	Status status(unionSample);
	getStatus(&status);

	std::vector<double> bet(status.k);
	std::vector<double> errL(status.k);
	std::vector<double> errU(status.k);

	computeBetErr(&status, bet, errL, errU);

	for (count i = 0; i < unionSample; ++i) {
		count v = status.top[i];
		approx[v] = approx[v] / (double)nPairs;
	}

	while (b - a > err / 10.) {
		c = (b + a) / 2.;
		sum = 0;
		for (count i = 0; i < unionSample; ++i) {
			sum += std::exp(-c * errL[i] * errL[i] / bet[i]);
			sum += std::exp(-c * errU[i] * errU[i] / bet[i]);
		}

		sum += std::exp(-c * errL[unionSample - 1] * errL[unionSample - 1] /
										bet[unionSample - 1]) *
					 (n - unionSample);
		sum += std::exp(-c * errU[unionSample - 1] * errU[unionSample - 1] /
										bet[unionSample - 1]) *
					 (n - unionSample);

		if (sum >= delta / 2. * (1 - balancingFactor)) {
			a = c;
		} else {
			b = c;
		}
	}

	delta_l_min_guess = std::exp(-b * errL[unionSample - 1] *
															 errL[unionSample - 1] / bet[unionSample - 1]) +
											delta * balancingFactor / 4. / (double)n;
	delta_u_min_guess = std::exp(-b * errU[unionSample - 1] *
															 errU[unionSample - 1] / bet[unionSample - 1]) +
											delta * balancingFactor / 4. / (double)n;

	std::fill(delta_l_guess.begin(), delta_l_guess.end(), delta_l_min_guess);
	std::fill(delta_u_guess.begin(), delta_u_guess.end(), delta_u_min_guess);

	for (count i = 0; i < unionSample; ++i) {
		node v = status.top[i];
		delta_l_guess[v] = std::exp(-b * errL[i] * errL[i] / bet[i]) +
											 delta * balancingFactor / 4. / (double)n;
	}
}

void KadabraBetweenness::init() {
	approx.resize(n);
	delta_l_guess.resize(n);
	delta_u_guess.resize(n);
	nPairs = 0;
}

void KadabraBetweenness::run() {
	Aux::Timer timer;
	timer.start();
	init();

	// TODO: check if exact diameter is required or if an upper bound is correct.
	// If so, which is the maximum relative error?
	Diameter diam(G, estimatedRange, 0);
	diam.run();
	// Getting diameter upper bound
	int32_t diameter = diam.getDiameter().second;

	omega =
			0.5 / err / err * (std::log2(diameter - 1) + 1 + std::log(0.5 / delta));

	const count tau = omega / startFactor;

	if (unionSample == 0) {
		unionSample =
				std::min(n, (count)std::max((2 * std::sqrt(G.numberOfEdges()) /
																		 omp_get_max_threads()),
																		k + 20.));
	}

	std::vector<count> seeds(omp_get_max_threads());
	for (count i = 0; i < omp_get_max_threads(); ++i) {
		seeds[i] = std::rand();
	}

#pragma omp parallel
	{
		SpSampler sampler(G, seeds[omp_get_thread_num()]);
		while (nPairs <= tau) {
			oneRound(sampler);
		}
	}

	timer.stop();
	std::cout << "Time for first part " << timer.elapsedMilliseconds() / 1000.
						<< std::endl;
	timer.start();
	computeDeltaGuess();
	timer.stop();
	std::cout << "Time for delta guess " << timer.elapsedMilliseconds() / 1000.
						<< std::endl;
	nPairs = 0;
	top.init();
	std::fill(approx.begin(), approx.end(), 0);
	double time_comp_finished = 0.;
	double time_round = 0.;
	double time_get_status = 0.;
#pragma omp parallel
	{
		SpSampler sampler(G, seeds[omp_get_thread_num()]);
		Status status(unionSample);
		status.nPairs = 0;
		bool stop = false;

		while (!stop && status.nPairs < omega) {
			timer.start();
			for (unsigned short i = 0; i <= 10; ++i) {
				oneRound(sampler);
			}
			timer.stop();
			time_round += timer.elapsedMilliseconds();
			timer.start();
			getStatus(&status);
			timer.stop();
			time_get_status += timer.elapsedMilliseconds();
			timer.start();
			stop = computeFinished(&status);
			timer.stop();
			time_comp_finished += timer.elapsedMilliseconds();
		}
	}

	nPairs += tau;
	Status status(unionSample);
	getStatus(&status);
	fillResult(&status);

	hasRun = true;
	std::cout << "Time for compute finished = " << time_comp_finished / 1000.
						<< std::endl;
	std::cout << "Time one round = " << time_round / 1000. << std::endl;
	std::cout << "Time get status = " << time_get_status / 1000. << std::endl;
}

SpSampler::SpSampler(const Graph &G, const count seed)
		: G(G), n(G.upperNodeIdBound()), pred(n, false, true) {
	q.resize(n);
	ballInd.resize(n);
	dist.resize(n);
	nPaths.resize(n);
}

std::vector<node> SpSampler::randomPath() {
	node u = G.randomNode();
	node v = G.randomNode();
	while (u == v) {
		v = G.randomNode();
	}

	count endQ = 2;
	q[0] = u;
	q[1] = v;

	ballInd[u] = 1;
	ballInd[v] = 2;

	dist[u] = 0;
	dist[v] = 0;
	nPaths[u] = 1;
	nPaths[v] = 1;

	std::vector<std::pair<node, node>> spEdges;

	node x, randomEdge;
	bool hasToStop = false;
	bool useDegreeIn;
	count startU = 0, startV = 1, endU = 1, endV = 2, startCur, endCur,
				*newEndCur;
	count sumDegsU = 0, sumDegsV = 0, *sumDegsCur;
	count totWeight = 0, curEdge = 0;

	auto procNeighbor = [&](node x, node y) {
		if (ballInd[y] == 0) {
			(*sumDegsCur) += getDegree(G, y, useDegreeIn);
			nPaths[y] = nPaths[x];
			assert(ballInd[x] != 0);
			ballInd[y] = ballInd[x];
			q[endQ++] = y;
			(*newEndCur)++;
			pred.addEdge(y, x);
			dist[y] = dist[x] + 1;
		} else if (ballInd[x] != ballInd[y]) {
			hasToStop = true;
			spEdges.push_back(std::make_pair(x, y));
		} else if (dist[y] == dist[x] + 1) {
			nPaths[y] += nPaths[x];
			pred.addEdge(y, x);
		}
	};

	while (!hasToStop) {
		if (sumDegsU <= sumDegsV) {
			startCur = startU;
			endCur = endU;
			startU = endQ;
			newEndCur = &endU;
			endU = endQ;
			sumDegsU = 0;
			sumDegsCur = &sumDegsU;
			useDegreeIn = false;
		} else {
			startCur = startV;
			endCur = endV;
			startV = endQ;
			newEndCur = &endV;
			endV = endQ;
			sumDegsV = 0;
			sumDegsCur = &sumDegsV;
			useDegreeIn = true;
		}

		while (startCur < endCur) {
			x = q[startCur++];

			if (useDegreeIn) {
				G.forInNeighborsOf(x, [&](node y) { procNeighbor(x, y); });
			} else {
				G.forNeighborsOf(x, [&](node y) { procNeighbor(x, y); });
			}
		}

		if (*sumDegsCur == 0) {
			hasToStop = true;
		}
	}

	if (spEdges.size() == 0) {
		removeAllEdges(endQ);
		assert(pred.numberOfEdges() == 0);
		std::fill(ballInd.begin(), ballInd.end(), 0);
		return std::vector<node>();
	}

	for (auto p : spEdges) {
		totWeight += nPaths[p.first] * nPaths[p.second];
	}

	randomEdge = Aux::Random::integer(totWeight);
	std::vector<node> path;

	for (auto p : spEdges) {
		curEdge += nPaths[p.first] * nPaths[p.second];
		if (curEdge > randomEdge) {
			backtrackPath(u, v, p.first, path);
			backtrackPath(u, v, p.second, path);
			break;
		}
	}

	std::fill(ballInd.begin(), ballInd.end(), 0);
	removeAllEdges(endQ);
	return path;
}

void SpSampler::backtrackPath(const node u, const node v, const node start,
															std::vector<node> &path) {
	if (start == u || start == v) {
		return;
	}

	count totWeight = nPaths[start];
	node randomPred, curPred = 0;
	node w = 0;

	path.push_back(start);
	randomPred = Aux::Random::integer(totWeight);
	assert((pred.neighbors(start)).size() > 0);

	for (node t : pred.neighbors(start)) {
		w = t;
		curPred += nPaths[v];
		if (curPred > randomPred) {
			break;
		}
	}

	if (w != u && w != v) {
		backtrackPath(u, v, w, path);
	}
}

void SpSampler::removeAllEdges(const count endQ) {
	std::vector<node> resizedQ(endQ);
	std::copy(q.begin(), q.begin() + endQ, resizedQ.begin());
	pred.removeEdgesFromIsolatedSet(resizedQ);
}

count SpSampler::getDegree(const Graph &graph, node z, bool useDegreeIn) {
	return useDegreeIn ? graph.degreeIn(z) : graph.degree(z);
}
} // namespace NetworKit
