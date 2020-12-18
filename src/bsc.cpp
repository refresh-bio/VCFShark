// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2020-12-18
// *******************************************************************************************

#include "bsc.h"

#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;

int CBSCWrapper::features;

// *******************************************************************************************
CBSCWrapper::CBSCWrapper()
{
}

// *******************************************************************************************
CBSCWrapper::~CBSCWrapper()
{
}

// *******************************************************************************************
bool CBSCWrapper::InitLibrary(int _features)
{
	features = _features;

	bsc_init(features);

	return true;
}

// *******************************************************************************************
bool CBSCWrapper::InitCompress(uint32_t _block_size, uint32_t _lzp_hash_size, uint32_t _lzp_min_len, uint32_t _coder)
{
	block_size = _block_size;
	lzp_hash_size = _lzp_hash_size;
	lzp_min_len = _lzp_min_len;
	coder = _coder;

	return true;
}

// *******************************************************************************************
bool CBSCWrapper::InitCompress(bsc_params_t params)
{
	block_size = params.block_size;
	lzp_hash_size = params.lzp_hash_size;
	lzp_min_len = params.lzp_min_len;
	coder = params.coder;

	return true;
}

// *******************************************************************************************
bool CBSCWrapper::InitDecompress()
{
	return true;
}

// *******************************************************************************************
bool CBSCWrapper::Compress(const vector<uint8_t>& v_input, vector<uint8_t>& v_output)
{
	const unsigned char* vi = (const unsigned char*)v_input.data();

	v_output.resize(v_input.size() + LIBBSC_HEADER_SIZE);

	unsigned char* vo = (unsigned char*)v_output.data();

	auto c_size = bsc_compress(vi, vo, (int) v_input.size(), lzp_hash_size, lzp_min_len, LIBBSC_BLOCKSORTER_BWT, coder, features);

	if (c_size == LIBBSC_NOT_COMPRESSIBLE)
		c_size = bsc_store(vi, vo, (int) v_input.size(), features);

	v_output.resize(c_size);
	v_output.shrink_to_fit();

	return true;
}

// *******************************************************************************************
bool CBSCWrapper::Decompress(vector<uint8_t>& v_input, vector<uint8_t>& v_output)
{
	vector<uint8_t> vi = v_input;
	vector<uint8_t> vo;

	int p_block_size;
	int p_data_size;

	const unsigned char* ci = (const unsigned char*)vi.data();
	bsc_block_info(ci, LIBBSC_HEADER_SIZE, &p_block_size, &p_data_size, 0);

#ifdef LOG_INFO
	cout << "Block size " << p_block_size << endl;
	cout << "Data size  " << p_data_size << endl;
#endif

	unsigned char* co = (unsigned char*)malloc(p_data_size);

	bsc_decompress(ci, p_block_size, co, p_data_size, 0);

	vo.assign(co, co + p_data_size);

	free(co);

	v_output = move(vo);

	return true;
}

// EOF
