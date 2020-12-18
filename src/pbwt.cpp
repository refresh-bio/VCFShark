// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2020-12-18
// *******************************************************************************************

#include "pbwt.h"
#include <numeric>
#include <algorithm>
#include <iterator>
#include "utils.h"

#include <iostream>

// ************************************************************************************
CPBWT::CPBWT()
{
}

// ************************************************************************************
CPBWT::~CPBWT()
{
}

// ************************************************************************************
void CPBWT::adjust_size(uint32_t new_size)
{
	size_t p_size = v_perm_prev.size();
	size_t c_size = new_size;

	if (c_size > p_size)
	{
		// If necessary to extend the permutation vector
		v_perm_prev.resize(c_size);

		if(v_removed_ids.size() < c_size)
			for (auto i = p_size; i < c_size; ++i)
				v_perm_prev[i] = (int)i;
		else
			for (auto i = p_size; i < c_size; ++i)
				v_perm_prev[i] = v_removed_ids[i];

	}
	else if (c_size < p_size)
	{
		// If necessary to compact the permutation vector
		if (v_removed_ids.size() < p_size)
			v_removed_ids.resize(p_size);

		copy_if(v_perm_prev.begin(), v_perm_prev.end(), v_removed_ids.begin() + c_size, [c_size](int x) {
			return x >= (int) c_size;
			});

		auto new_end = remove_if(v_perm_prev.begin(), v_perm_prev.end(), [c_size](int x) {
			return x >= (int) c_size;
			});

		v_perm_prev.erase(new_end, v_perm_prev.end());

	}
}

// ************************************************************************************
bool CPBWT::StartForward(const size_t _no_items, const size_t _neglect_limit)
{
	no_items = _no_items;
	neglect_limit = _neglect_limit;

	v_perm_cur.resize(no_items);
	v_tmp.resize(no_items);

	iota(v_perm_cur.begin(), v_perm_cur.end(), 0);
	v_perm_prev = v_perm_cur;
	
	return true;
}

// ************************************************************************************
bool CPBWT::StartReverse(const size_t _no_items, const size_t _neglect_limit)
{
	no_items = _no_items;
	neglect_limit = _neglect_limit;

	v_perm_cur.resize(no_items);

	iota(v_perm_cur.begin(), v_perm_cur.end(), 0);
	v_perm_prev = v_perm_cur;
	
	v_tmp.clear();
	v_tmp.resize(no_items, 0u);

	return true;
}

// ************************************************************************************
// Forward PBWT for non-binary alphabet
bool CPBWT::EncodeFlexible(const uint32_t max_val, vector<uint32_t> &v_input, vector<pair<uint32_t, uint32_t>> &v_rle)
{
	
	v_hist.resize(max_val + 1);
	size_t c_size = v_input.size();

	uint32_t max_count;

	// Determine histogram of symbols
	calc_cumulate_histogram(v_input, v_hist, max_count);

	vector<int> v_perm_prev0;

	if (c_size != v_perm_prev.size())
	{
		if (c_size - max_count < neglect_limit)
			v_perm_prev0 = v_perm_prev;

		adjust_size((uint32_t) c_size);
	}
	else
		v_perm_prev0.clear();

	v_perm_cur.resize(c_size);

	uint8_t prev_symbol = (uint8_t) v_input[v_perm_prev[0]];
	uint32_t run_len = 0;

	v_rle.clear();

	// Make PBWT
	for (size_t i = 0; i < c_size; ++i)
	{
		uint8_t cur_symbol = (uint8_t) v_input[v_perm_prev[i]];

		if (cur_symbol == prev_symbol)
			++run_len;
		else
		{
			v_rle.emplace_back(prev_symbol, run_len);
			prev_symbol = cur_symbol;
			run_len = 1;
		}

		v_perm_cur[v_hist[cur_symbol]] = v_perm_prev[i];
		++v_hist[cur_symbol];
	}

	v_rle.emplace_back(prev_symbol, run_len);

	// Swap only if no. of non-zeros is larger than neglect_limit
	if (c_size - max_count >= neglect_limit)
		swap(v_perm_prev, v_perm_cur);
	else if(!v_perm_prev0.empty())
		v_perm_prev = v_perm_prev0;

	return true;
}

// ************************************************************************************
// Reverse PBWT for non-binary alphabet
bool CPBWT::DecodeFlexible(const uint32_t max_val, const vector<pair<uint32_t, uint32_t>>& v_rle, vector<uint32_t>& v_output)
{
	vector<uint32_t> v_hist(max_val + 1);
	uint32_t max_count;

	uint32_t no_items = 0;

	for (auto x : v_rle)
		no_items += x.second;

	v_output.resize(no_items);

	calc_cumulate_histogram(v_rle, v_hist, max_count);

	size_t c_size = no_items;

	vector<int> v_perm_prev0;

	if (c_size != v_perm_prev.size())
	{
		if (c_size - max_count < neglect_limit)
			v_perm_prev0 = v_perm_prev;

		adjust_size((uint32_t) c_size);
	}
	else
		v_perm_prev0.clear();

	auto p_rle = v_rle.begin();
	uint8_t cur_symbol = (uint8_t) p_rle->first;
	uint32_t cur_cnt = p_rle->second;

	v_perm_cur.resize(no_items);

	// Make PBWT
	for (size_t i = 0; i < no_items; ++i)
	{
		v_output[v_perm_prev[i]] = cur_symbol;

		v_perm_cur[v_hist[cur_symbol]] = v_perm_prev[i];
		++v_hist[cur_symbol];

		if (--cur_cnt == 0)
		{
			++p_rle;
			if (i + 1 < no_items)
			{
				cur_symbol = (uint8_t) p_rle->first;
				cur_cnt = p_rle->second;
			}
		}
	}

	// Swap only if no. of non-zeros is larger than neglect_limit
	if (no_items - max_count >= neglect_limit)
		swap(v_perm_prev, v_perm_cur);
	else if (!v_perm_prev0.empty())
		v_perm_prev = v_perm_prev0;

	return true;
}

// EOF
