// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include "buffer.h"

#include <iostream>

// ************************************************************************************
CBuffer::CBuffer()
{
	type = buffer_t::none;
	max_size = 0;

	v_size_pos = 0;
	v_data_pos = 0;

	is_function = false;
	is_no_data = false;
}

// ************************************************************************************
CBuffer::~CBuffer()
{
}

// ************************************************************************************
void CBuffer::SetMaxSize(uint32_t _max_size, uint32_t _offset)
{
	max_size = _max_size;
#ifndef IGNORE_OFFSET_IN_FIRST_BLOCK
	offset = _offset;
#endif
}

// ************************************************************************************
void CBuffer::GetBuffer(vector<uint32_t>& _v_size, vector<uint8_t>& _v_data)
{
	v_size.shrink_to_fit();
	v_data.shrink_to_fit();

	swap(_v_size, v_size);
	swap(_v_data, v_data);

	v_size.clear();
	v_data.clear();

/*	v_size.shrink_to_fit();
	v_data.shrink_to_fit();

	v_size.reserve(_v_size.size() / 4);
	v_data.reserve(_v_data.size() / 4);*/

	type = buffer_t::none;

#ifndef IGNORE_OFFSET_IN_FIRST_BLOCK
	offset = 0;
#endif
}

// ************************************************************************************
void CBuffer::permute_integer_series_forward()
{
	int no_items = 0;
	int no_series = 0;
	int series_size = 0;

	for (auto x : v_size)
	{
		if (x)
		{
			if (series_size == 0)
				series_size = x;
			else
				series_size = gcd(x, series_size);

			no_items += x;
		}
	}

	if (series_size <= 1)
	{
		v_data.insert(v_data.begin(), 0);
		return;
	}

	no_series = no_items / series_size;

	vector<uint32_t> v_tmp(no_series * series_size);

	auto v_data0 = v_data;		// Copy of v_data

	uint8_t* p = v_data.data();
	uint32_t val;

	uint32_t p_val = 0;
	int no_orig = 0;
	int no_perm = 0;

	for (int i = 0; i < no_series; ++i)
		for(int j = 0; j < series_size; ++j)
		{
			p += decode_var_int(p, val);
			v_tmp[j * no_series + i] = val;

			if (p_val != val)
			{
				++no_orig;
				p_val = val;
			}
		}

	v_data.clear();
	
	char* q = (char*) v_tmp.data();
	p_val = 0;

	v_data.emplace_back(1);		// Marker of permutation

	for (int i = 0; i < no_series * series_size; ++i)
	{
		encode_var_int(q + 4 * i);

		uint32_t val = v_tmp[i];

		if (p_val != val)
		{
			++no_perm;
			p_val = val;
		}
	}

	if (no_orig < no_perm)
	{
		v_data = move(v_data0);
		v_data.insert(v_data.begin(), 0);
	}
}

// ************************************************************************************
void CBuffer::permute_integer_series_backward()
{
	if (v_data.front() == 0)			// no permutation
	{
		v_data.erase(v_data.begin());
		
		return;
	}

	int no_items = 0;
	int no_series = 0;
	int series_size = 0;

	for (auto x : v_size)
	{
		if (x)
		{
			if (series_size == 0)
				series_size = x;
			else
				series_size = gcd(x, series_size);

			no_items += x;
		}
	}
		
	no_series = no_items / series_size;

	vector<uint32_t> v_tmp(no_series * series_size);

	uint8_t* p = v_data.data() + 1;
	uint32_t val;

	for (int i = 0; i < series_size; ++i)
		for (int j = 0; j < no_series; ++j)
		{
			p += decode_var_int(p, val);
			v_tmp[j * series_size + i] = val;
		}

	v_data.clear();

	char* q = (char*)v_tmp.data();

	for (int i = 0; i < no_series * series_size; ++i)
		encode_var_int(q + 4 * i);
}

// ************************************************************************************
void CBuffer::permute_float_series_forward()
{
	int no_items = 0;
	int no_series = 0;
	int series_size = 0;

	for (auto x : v_size)
	{
		if (x)
		{
			if (series_size == 0)
				series_size = x;
			else
				series_size = gcd(x, series_size);

			no_items += x;
		}
	}

	if (series_size <= 1)
	{
		v_data.insert(v_data.begin(), 0);
		return;
	}

	no_series = no_items / series_size;

	vector<float> v_tmp(no_series * series_size);

	auto v_data0 = v_data;		// Copy of v_data

	char* p_src = (char*) v_data0.data();
	char* p_dest = (char*) v_data.data();

	for (int i = 0; i < no_series; ++i)
		for (int j = 0; j < series_size / 3; ++j)
			memcpy(p_dest + j * no_series * 12 + i * 12, p_src + i * series_size * 4 + j * 12, 12);

	return;

}

// ************************************************************************************
void CBuffer::permute_float_series_backward()
{
	if (v_data.front() == 0)			// no permutation
	{
		v_data.erase(v_data.begin());
		
		return;
	}

	int no_items = 0;
	int no_series = 0;
	int series_size = 0;

	for (auto x : v_size)
	{
		if (x)
		{
			if (series_size == 0)
				series_size = x;
			else
				series_size = gcd(x, series_size);

			no_items += x;
		}
	}
		
	no_series = no_items / series_size;

	vector<uint32_t> v_tmp(no_series * series_size);

	uint8_t* p = v_data.data() + 1;
	uint32_t val;

	for (int i = 0; i < series_size; ++i)
		for (int j = 0; j < no_series; ++j)
		{
			p += decode_var_int(p, val);
			v_tmp[j * series_size + i] = val;
		}

	v_data.clear();

	char* q = (char*)v_tmp.data();

	for (int i = 0; i < no_series * series_size; ++i)
		encode_var_int(q + 4 * i);
}

// ************************************************************************************
void CBuffer::FuncInt(char*& p, uint32_t& size, char* src_p, uint32_t src_size)
{
	if (!src_size)
	{
		size = 0;
		p = nullptr;
		
		return;
	}

	if (fun.empty())		// identity
	{
		p = new char[src_size * 4];
		copy_n(src_p, src_size * 4, p);
		size = src_size;
	}
	else
	{
		vector<uint8_t> src_vec(src_p, src_p + src_size);

		auto& dest_vec = fun[src_vec];

		p = new char[dest_vec.size()];
		copy_n(dest_vec.data(), dest_vec.size(), p);
		size = (uint32_t) dest_vec.size();
	}
}

// ************************************************************************************
void CBuffer::FuncReal(char*& p, uint32_t& size, char* src_p, uint32_t src_size)
{
	if (!src_size)
	{
		size = 0;
		p = nullptr;

		return;
	}

	if (fun.empty())	// identity
	{
		p = new char[src_size * 4];
		copy_n(src_p, src_size * 4, p);
		size = src_size;
	}
	else
	{
		vector<uint8_t> src_vec(src_p, src_p + src_size);

		auto& dest_vec = fun[src_vec];

		p = new char[dest_vec.size()];
		copy_n(dest_vec.data(), dest_vec.size(), p);
		size = (uint32_t) dest_vec.size();
	}
}

// ************************************************************************************
void CBuffer::SetBuffer(vector<uint32_t>& _v_size, vector<uint8_t>& _v_data)
{
	v_size = move(_v_size);
	v_data = move(_v_data);

	_v_size.clear();
	_v_data.clear();

	v_size_pos = 0;
	v_data_pos = 0;

	is_no_data = v_size.empty();
}

// ************************************************************************************
void CBuffer::SetFunction(function_data_item_t& _fun)
{
	fun = _fun;

	is_function = true;
}

// EOF
