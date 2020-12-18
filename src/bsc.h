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
#include "defs.h"

#include <libbsc.h>

using namespace std;

// ************************************************************************************
struct bsc_params_t {
	uint32_t block_size;
	uint32_t lzp_hash_size;
	uint32_t lzp_min_len;
	uint32_t coder;
};

// ************************************************************************************
class CBSCWrapper
{
	uint32_t block_size;
	uint32_t lzp_hash_size;
	uint32_t lzp_min_len;
	uint32_t coder;

	static int features;

public:
	CBSCWrapper();
	~CBSCWrapper();

	static bool InitLibrary(int _features);

	bool InitDecompress();
	bool InitCompress(uint32_t _block_size, uint32_t _lzp_hash_size, uint32_t _lzp_min_len, uint32_t _coder);
	bool InitCompress(bsc_params_t params);

	bool Compress(const vector<uint8_t>& v_input, vector<uint8_t>& v_output);
	static bool Decompress(vector<uint8_t>& v_input, vector<uint8_t>& v_output);
};

// EOF
