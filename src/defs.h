#pragma once
// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2020-12-18
// *******************************************************************************************

#include <array>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <vector>

#define OUR_STRTOL


using namespace std;

#ifndef uint32_t
typedef unsigned int uint32_t;
#endif



#ifndef uint8_t
typedef unsigned char uint8_t;
#endif


typedef array<uint32_t, 4> stats_t;
typedef array<uint64_t, 4> stats64_t;

typedef uint64_t context_t;

typedef array<pair<uint8_t, uint32_t>, 2> run_t;

const uint32_t SIGMA = 4u;

template <typename T>
class pair_hash
{
public:
	size_t operator()(const pair<T, T>& x) const {
		return hash<T>{}(x.first) ^ hash<T>{}(x.second);
	}
};

template <typename T>
class vector_hash
{
public:
	size_t operator()(const vector<T>& x) const {
		size_t r = 0;
		for (auto& c : x)
			r ^= hash<T>{}(c);

		return r;
	}
};

using function_size_item_t = map<uint32_t, uint32_t>;
using function_size_graph_t = map<pair<int, int>, function_size_item_t>;

using function_data_item_t = unordered_map<vector<uint8_t>, vector<uint8_t>, vector_hash<uint8_t>>;
using function_data_graph_t = unordered_map<pair<int, int>, function_data_item_t, pair_hash<int>>;

// EOF
