/*
 * PQVector.h
 *
 *  Created on: 21.09.2018
 *      Author: Eugenio Angriman
 */

#ifndef PQVECTOR_H_
#define PQVECTOR_H_

#include <algorithm>
#include <vector>

#include "../auxiliary/Log.h"

namespace Aux {
/*
 * A humble attempt to have a faster prio queue for Kadabra.
 */
class PQVector {
private:
	std::vector<uint64_t> elements;
	std::vector<double> values;
	std::vector<uint64_t> position;
	uint64_t virtual_size;
	const uint64_t size;
	const uint64_t max_key;

public:
	PQVector(const uint64_t size, const uint64_t max_key);
	void insert(const uint64_t new_element, const double new_value);
	double get_value(const uint64_t i) const { return values[i]; }
	uint64_t get_element(const uint64_t i) const { return elements[i]; }
	uint64_t get_size() const { return virtual_size; }
	void print() const {
		INFO("Elements = ", elements);
		INFO("Values   = ", values);
		INFO("Position = ", position);
		INFO("Size = ", virtual_size);
	}
	void init();
};

inline Aux::PQVector::PQVector(const uint64_t size, const uint64_t max_key)
		: size(size), max_key(max_key) {
	init();
}

inline void Aux::PQVector::init() {
	elements.resize(size);
	values.assign(size, 0.);
	position.resize(max_key);
	for (uint64_t i = 0; i < size; ++i) {
		elements[i] = i;
		position[i] = i;
	}
	virtual_size = 0;
}

inline void Aux::PQVector::insert(const uint64_t new_element,
																	const double new_value) {
	uint64_t ub = std::upper_bound(values.begin(), values.begin() + virtual_size,
																 new_value, std::greater<double>()) -
								values.begin();
	uint64_t old_pos;
	// We assume that if the same key is inserted again, its value will be
	// greater than before (i.e., old_pos > ub).
	if (ub < size) {
		old_pos = position[new_element];
		if (virtual_size < size) {
			++virtual_size;
		}
	} else {
		return;
	}
	++old_pos;
	std::rotate(values.begin() + ub, values.begin() + old_pos - 1,
							values.begin() + old_pos);
	std::rotate(elements.begin() + ub, elements.begin() + old_pos - 1,
							elements.begin() + old_pos);
	values[ub] = new_value;
	elements[ub] = new_element;
	position[new_element] = ub;
	for (auto it = elements.begin() + ub + 1; it < elements.begin() + old_pos;
			 ++it) {
		++position[*it];
	}
}

} // namespace Aux

#endif
