#pragma once
// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2020-12-18
// *******************************************************************************************

#include <vector>
#include <set>
#include <map>
#include <utility>
#include <cstdint>

using namespace std;

// ************************************************************************************
class CGraphOptimizer
{
public:
	struct node_t {
		int id;
		uint64_t cost;
		node_t(int _id, uint64_t _cost) : id(_id), cost(_cost) {};
	};

	struct edge_t {
		int id_from;
		int id_to;
		bool equality;
		uint64_t cost;
		edge_t(int _id_from, int _id_to, bool _equality, uint64_t _cost) : id_from(_id_from), id_to(_id_to), equality(_equality), cost(_cost) {};
	};

	struct edge_id_t {
		int from;
		int to;
		edge_id_t() : from(-1), to(-1) {};
		edge_id_t(int _from, int _to) : from(_from), to(_to) {};
	};

private:

	struct node_data_t {
		uint64_t cost;
		int in_degree;
		int out_degree;
		node_data_t() : cost(0), in_degree(0), out_degree(0) {};
		node_data_t(uint64_t _cost, int _in_degree, int _out_degree) : cost(_cost), in_degree(_in_degree), out_degree(_out_degree) {};
	};

	struct edge_data_t {
		bool equality;
		uint64_t cost;
		edge_data_t() : equality(false), cost(0) {};
		edge_data_t(bool _equality, uint64_t _cost) : equality(_equality), cost(_cost) {};
	};

	// Graph
	map<int, node_data_t> m_nodes;		
	map<edge_id_t, edge_data_t> m_edges;

	// Nodes and edges for output
	vector<int> v_nodes;
	vector<pair<int, int>> v_edges;

	void load_graph(const vector<node_t>& in_nodes, const vector<edge_t>& in_edges);
	bool remove_isolated_nodes();
	bool remove_equality_edges();
	bool remove_expensive_edges();

	bool remove_best_edge();
	bool remove_random_edge();

	bool strip_edges();

	void remove_edges_with_target(int to);

public:
	CGraphOptimizer() {};
	~CGraphOptimizer() {};

	bool Optimize(const vector<node_t>& in_nodes, const vector<edge_t>& in_edges, vector<pair<int, bool>> &out_nodes, vector<pair<int, int>> &out_edges);
};

bool operator<(const CGraphOptimizer::edge_id_t x, const CGraphOptimizer::edge_id_t y);

// EOF
