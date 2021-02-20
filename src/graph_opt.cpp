// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include "graph_opt.h"
#include <unordered_set>
#include <cassert>
#include <limits>

// ******************************************************************************
bool CGraphOptimizer::Optimize(const vector<node_t>& in_nodes, const vector<edge_t>& in_edges, vector<pair<int, bool>> &out_nodes, vector<pair<int, int>> &out_edges)
{
	load_graph(in_nodes, in_edges);
	remove_isolated_nodes();
	remove_equality_edges();
	remove_isolated_nodes();

	m_edges.clear();
	for (auto& n : m_nodes)
		n.second.in_degree = n.second.out_degree = 0;

	remove_expensive_edges();
	remove_isolated_nodes();
	strip_edges();
	remove_isolated_nodes();

	out_nodes.clear();

	for (auto p = v_nodes.rbegin(); p != v_nodes.rend(); ++p)
	{
		bool explicite = true;

		for(auto &e : v_edges)
			if (e.second == *p)
			{
				explicite = false;
				break;
			}

		out_nodes.emplace_back(*p, explicite);
	}
	
	out_edges = move(v_edges);

	return true;
}

// ******************************************************************************
void CGraphOptimizer::load_graph(const vector<node_t>& in_nodes, const vector<edge_t>& in_edges)
{
	m_nodes.clear();
	m_edges.clear();

	v_nodes.clear();
	v_nodes.clear();

	for (auto& n : in_nodes)
		m_nodes[n.id] = node_data_t(n.cost, 0, 0);

	for (auto& e : in_edges)
	{
		m_edges[edge_id_t(e.id_from, e.id_to)] = edge_data_t(e.equality, e.cost);
		m_nodes[e.id_from].out_degree++;
		m_nodes[e.id_to].in_degree++;
	}
}

// ******************************************************************************
bool CGraphOptimizer::remove_isolated_nodes()
{
	vector<int> nodes_to_remove;

	for (auto& n : m_nodes)
		if (n.second.in_degree == 0 && n.second.out_degree == 0)
			nodes_to_remove.emplace_back(n.first);

	for (auto& n_id : nodes_to_remove)
	{
		v_nodes.emplace_back(n_id);
		m_nodes.erase(n_id);
	}

	return !nodes_to_remove.empty();
}

// ******************************************************************************
bool CGraphOptimizer::remove_equality_edges()
{
	unordered_set<int> nodes_to_remove;
	vector<edge_id_t> edges_to_remove;

	// Look for equality edges and mark nodes to remove
	for (auto& n1 : m_nodes)
		if (nodes_to_remove.count(n1.first))
			continue;
		else
			for(auto &n2 : m_nodes)
				if (nodes_to_remove.count(n2.first))
					continue;
				else
				{
					auto p = m_edges.find(edge_id_t(n1.first, n2.first));
					if (p != m_edges.end() && p->second.equality)
					{
						nodes_to_remove.insert(n2.first);
						v_edges.emplace_back(make_pair(n1.first, n2.first));
						v_nodes.emplace_back(n2.first);
					}
				}

	// Mark edges for nodes_to_remove to remove
	for (auto& e : m_edges)
		if (nodes_to_remove.count(e.first.from) || nodes_to_remove.count(e.first.to))
			edges_to_remove.emplace_back(e.first);

	// Remove redundant edges
	for (auto& e_id : edges_to_remove)
	{
		m_edges.erase(e_id);
		m_nodes[e_id.from].out_degree--;
		m_nodes[e_id.to].in_degree--;
	}

	// Remove redundant nodes
	for (auto& n_id : nodes_to_remove)
		m_nodes.erase(n_id);

	return !nodes_to_remove.empty();
}

// ******************************************************************************
// Remove edges of to high cost
bool CGraphOptimizer::remove_expensive_edges()
{
	vector<edge_id_t> edges_to_remove;
	vector<int> nodes_to_remove;

	for (auto& e : m_edges)
		if (e.second.cost > m_nodes[e.first.to].cost)
			edges_to_remove.emplace_back(e.first);

	for (auto& e : edges_to_remove)
	{
		m_edges.erase(e);
		m_nodes[e.from].out_degree--;
		m_nodes[e.to].in_degree--;
	}

	return !edges_to_remove.empty();
}

// ******************************************************************************
// Remove edge with the largest benefit
bool CGraphOptimizer::remove_best_edge()
{
	edge_id_t best_edge{ -1, -1 };
	int64_t best_gain = 0;

	for(auto &e : m_edges)
		if (m_nodes[e.first.to].out_degree == 0)
		{
			int64_t cur_gain = (int64_t)m_nodes[e.first.to].cost - (int64_t)e.second.cost;

			if (cur_gain > best_gain)
			{
				best_gain = cur_gain;
				best_edge = e.first;
			}
		}

	if (best_gain > 0)	// If anythink found
	{
		v_nodes.emplace_back(best_edge.to);
		v_edges.emplace_back(best_edge.from, best_edge.to);

		remove_edges_with_target(best_edge.to);

		assert(m_nodes[best_edge.to].in_degree == 0);
		m_nodes.erase(best_edge.to);

		return true;
	}
	else
		return false;
}

// ******************************************************************************
// Remove ''random'' edge - helpful; for cyclic relationships
bool CGraphOptimizer::remove_random_edge()
{
	edge_id_t worse_edge{ -1, -1 };
	int64_t worse_gain = numeric_limits<int64_t>::max();

	for (auto& e : m_edges)
	{
		int64_t cur_gain = (int64_t)m_nodes[e.first.to].cost - (int64_t)e.second.cost;

		if (cur_gain < worse_gain)
		{
			worse_gain = cur_gain;
			worse_edge = e.first;
		}
	}

	assert(worse_edge.from != -1);

	m_edges.erase(worse_edge);
	m_nodes[worse_edge.from].out_degree--;
	m_nodes[worse_edge.to].in_degree--;

	return true;
}

// ******************************************************************************
bool CGraphOptimizer::strip_edges()
{
	while (!m_edges.empty())
		if (!remove_best_edge())
			remove_random_edge();

	return true;
}

// ******************************************************************************
void CGraphOptimizer::remove_edges_with_target(int to)
{
	vector<edge_id_t> edges_to_remove;

	for (auto& e : m_edges)
		if (e.first.to == to)
			edges_to_remove.emplace_back(e.first);

	for (auto& e_id : edges_to_remove)
	{
		m_edges.erase(e_id);
		m_nodes[e_id.from].out_degree--;
		m_nodes[e_id.to].in_degree--;
	}
}

// ******************************************************************************
bool operator<(const CGraphOptimizer::edge_id_t x, const CGraphOptimizer::edge_id_t y)
{
	return x.from != y.from ? x.from < y.from : x.to < y.to;
}

// EOF
