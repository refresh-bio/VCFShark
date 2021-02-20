#pragma once
// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include <vector>
#include "defs.h"

using namespace std;

// ************************************************************************************
class CPBWT
{
	size_t no_items;
	size_t neglect_limit;

	vector<int> v_perm_cur;
	vector<int> v_perm_prev;
	vector<uint8_t> v_tmp;

	vector<int> v_removed_ids;

	vector<uint32_t> v_hist;
	vector<uint32_t> v_hist_complete;

	void adjust_size(uint32_t new_size);

public:
	CPBWT();
	~CPBWT();

	bool StartForward(const size_t _no_items, const size_t _neglect_limit);
	bool StartReverse(const size_t _no_items, const size_t _neglect_limit);

	bool EncodeFlexible(const uint32_t max_val, vector<uint32_t> &v_input, vector<pair<uint32_t, uint32_t>> &v_rle);
	bool DecodeFlexible(const uint32_t max_val, const vector<pair<uint32_t, uint32_t>> &v_rle, vector<uint32_t> &v_output);
};

// EOF
