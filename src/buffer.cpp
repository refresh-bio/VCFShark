// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2020-12-18
// *******************************************************************************************

#include "buffer.h"

#include <algorithm>
#include <cstring>
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
void CBuffer::SetMaxSize(uint32_t _max_size)
{
	max_size = _max_size;
}

// ************************************************************************************
bool CBuffer::IsFull(void)
{
	return v_data.size() + 4 * v_size.size() >= max_size;
}

// ************************************************************************************
void CBuffer::WriteFlag(uint8_t flag)
{
	v_size.emplace_back(flag);

	type = buffer_t::flag;
}

// ************************************************************************************
void CBuffer::WriteInt(char* p, uint32_t size)
{
	v_size.emplace_back(size);

	v_data.insert(v_data.end(), p, p + 4 * size);

	type = buffer_t::integer;
}

// ************************************************************************************
void CBuffer::WriteIntVarSize(char* p, uint32_t size)
{
	v_size.emplace_back(size);

	for (size_t i = 0; i < 4u * size; i += 4u)
		encode_var_int(p + i);
}

// ************************************************************************************
uint32_t CBuffer::encode_var_int(char* p)
{
	uint8_t* q = (uint8_t*)p;

	uint32_t val = 0;

	val += q[3];		val <<= 8;
	val += q[2];		val <<= 8;
	val += q[1];		val <<= 8;
	val += q[0];

	int32_t i_val = (int32_t)val;

	if (val == 0)
		v_data.emplace_back(0);
	else if (val == 0x80000000u)
		v_data.emplace_back(1);
	else if (i_val > 0 && i_val < 125)
		v_data.emplace_back(i_val + 1);
	else if (i_val < 0 && i_val > -125)
		v_data.emplace_back(i_val + 250);
	else if (i_val > 0 && i_val < 256 * 256)
	{
		v_data.emplace_back(250);
		v_data.emplace_back(i_val >> 8);
		v_data.emplace_back(i_val & 0xff);
	}
	else if (-i_val > 0 && -i_val < 256 * 256)
	{
		i_val = -i_val;
		v_data.emplace_back(251);
		v_data.emplace_back(i_val >> 8);
		v_data.emplace_back(i_val & 0xff);
	}
	else if (i_val > 0 && i_val < 256 * 256 * 256)
	{
		v_data.emplace_back(252);
		v_data.emplace_back(i_val >> 16);
		v_data.emplace_back((i_val >> 8) & 0xff);
		v_data.emplace_back(i_val & 0xff);
	}
	else if (-i_val > 0 && -i_val < 256 * 256 * 256)
	{
		i_val = -i_val;
		v_data.emplace_back(253);
		v_data.emplace_back(i_val >> 16);
		v_data.emplace_back((i_val >> 8) & 0xff);
		v_data.emplace_back(i_val & 0xff);
	}
	else if (i_val > 0)
	{
		v_data.emplace_back(254);
		v_data.emplace_back(i_val >> 24);
		v_data.emplace_back((i_val >> 16) & 0xff);
		v_data.emplace_back((i_val >> 8) & 0xff);
		v_data.emplace_back(i_val & 0xff);
	}
	else
	{
		i_val = -i_val;
		v_data.emplace_back(255);
		v_data.emplace_back(i_val >> 24);
		v_data.emplace_back((i_val >> 16) & 0xff);
		v_data.emplace_back((i_val >> 8) & 0xff);
		v_data.emplace_back(i_val & 0xff);
	}

	return 0;
}

// ************************************************************************************
uint32_t CBuffer::decode_var_int(uint8_t* p, uint32_t &val)
{
	int code = p[0];

	val = 0;

	if (code == 0)
		return 1;
	else if (code == 1)
	{
		val = 0x80000000u;
		return 1;
	}
	else if (code < 126)
	{
		val = code - 1;
		return 1;
	}
	else if (code < 250)
	{
		val = (uint32_t) (((int)code) - 250);
		return 1;
	}
	else if (code == 250)
	{
		val += (uint32_t)p[1];		val <<= 8;
		val += (uint32_t)p[2];
		return 3;
	}
	else if (code == 251)
	{
		val += (uint32_t)p[1];		val <<= 8;
		val += (uint32_t)p[2];
		
		val = (uint32_t) -((int)val);

		return 3;
	}
	else if (code == 252)
	{
		val += (uint32_t)p[1];		val <<= 8;
		val += (uint32_t)p[2];		val <<= 8;
		val += (uint32_t)p[3];
		return 4;
	}
	else if (code == 253)
	{
		val += (uint32_t)p[1];		val <<= 8;
		val += (uint32_t)p[2];		val <<= 8;
		val += (uint32_t)p[3];
		
		val = (uint32_t) -((int)val);

		return 4;
	}
	else if (code == 254)
	{
		val += (uint32_t)p[1];		val <<= 8;
		val += (uint32_t)p[2];		val <<= 8;
		val += (uint32_t)p[3];		val <<= 8;
		val += (uint32_t)p[4];
		return 5;
	}
	else if (code == 255)
	{
		val += (uint32_t)p[1];		val <<= 8;
		val += (uint32_t)p[2];		val <<= 8;
		val += (uint32_t)p[3];		val <<= 8;
		val += (uint32_t)p[4];
		
		val = (uint32_t) -((int)val);

		return 5;
	}

	return 0;	// !!! Never should be here
}

// ************************************************************************************
uint32_t CBuffer::encode_float(char* p)
{
	v_data.emplace_back(p[0]);
	v_data.emplace_back(p[1]);
	v_data.emplace_back(p[2]);
	v_data.emplace_back(p[3]);

	return 4;
}

// ************************************************************************************
uint32_t CBuffer::decode_float(uint8_t* p, float& val)
{
	uint32_t x = 0;

	x += p[3];		x <<= 8;
	x += p[2];		x <<= 8;
	x += p[1];		x <<= 8;
	x += p[0];

	memcpy(&val, &x, 4);

	return 4;
}

// ************************************************************************************
void CBuffer::WriteReal(char* p, uint32_t size)
{
	v_size.emplace_back(size);

	v_data.insert(v_data.end(), p, p + 4 * size);

	type = buffer_t::real;
}

// ************************************************************************************
void CBuffer::WriteText(char* p, uint32_t size)
{
	v_size.push_back(size);
	
	if(size)
		v_data.insert(v_data.end(), p, p + size);

	type = buffer_t::text;
}

// ************************************************************************************
void CBuffer::WriteInt64(int64_t x)
{
	int sign = 0;

	if (x < 0)
	{
		sign = 1;
		x = -x;
	}

	uint8_t bytes[8];
	int no_bytes = 0;
	size_t tmp = (size_t)x;

	for (; tmp; ++no_bytes)
	{
		bytes[no_bytes] = tmp & 0xff;
		tmp >>= 8;
	}

	v_size.emplace_back(sign + no_bytes * 2);

	for (int i = no_bytes - 1; i >= 0; --i)
		v_data.emplace_back(bytes[i]);
}

// ************************************************************************************
void CBuffer::ReadInt64(int64_t& x)
{
	int sign = v_size[v_size_pos++];

	int no_bytes = sign / 2;
	sign &= 1;

	x = 0;

	for (int i = 0; i < no_bytes; ++i)
		x += ((int64_t)v_data[v_data_pos++]) << (8 * (no_bytes - i - 1));

	if (sign)
		x = -x;
}

// ************************************************************************************
void CBuffer::GetBuffer(vector<uint32_t>& _v_size, vector<uint8_t>& _v_data)
{
	_v_size = move(v_size);
	_v_data = move(v_data);

	v_size.clear();
	v_data.clear();

	v_size.shrink_to_fit();
	v_data.shrink_to_fit();

	type = buffer_t::none;
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
		{
			memcpy(p_dest + j * no_series * 12 + i * 12, p_src + i * series_size * 4 + j * 12, 12);

		}

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
void CBuffer::ReadFlag(uint8_t &flag)
{
	if (v_size.empty())
	{
		flag = 0;

		return;
	}

	flag = (uint8_t) v_size[v_size_pos++];
}

// ************************************************************************************
void CBuffer::ReadInt(char*& p, uint32_t& size)
{
	if (v_size.empty())
	{
		p = nullptr;

		return;
	}

	size = v_size[v_size_pos++];

	if (size)
	{
		p = new char[size * 4];

		copy_n(v_data.begin() + v_data_pos, 4 * size, p);
		v_data_pos += 4 * size;
	}
	else
		p = nullptr;
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
void CBuffer::ReadIntVarSize(char*& p, uint32_t& size)
{
	if (v_size.empty())
	{
		p = nullptr;

		return;
	}

	size = v_size[v_size_pos++];

	if (size)
	{
		p = new char[size * 4];
		uint32_t val;

		uint8_t* q = v_data.data();

		for (size_t i = 0; i < size; ++i)
		{
			v_data_pos += decode_var_int(q + v_data_pos, val);
			memcpy(p + 4 * i, &val, 4);
		}
	}
	else
		p = nullptr;
}

// ************************************************************************************
void CBuffer::ReadReal(char*& p, uint32_t& size)
{
	if (v_size.empty())
	{
		p = nullptr;

		return;
	}

	size = v_size[v_size_pos++];

	if (size)
	{
		p = new char[size * 4];

		copy_n(v_data.begin() + v_data_pos, 4 * size, p);
		v_data_pos += 4 * size;
	}
	else
		p = nullptr;
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
void CBuffer::ReadText(char*& p, uint32_t& size)
{
	if (v_size.empty())
	{
		p = nullptr;

		return;
	}

	size = v_size[v_size_pos++];

	if (size)
	{
		p = new char[size + 1];

		copy_n(v_data.begin() + v_data_pos, size, p);
		p[size] = 0;
		v_data_pos += size;
	}
	else
		p = nullptr;
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

// ************************************************************************************
bool CBuffer::IsEmpty()
{
	return v_size_pos >= v_size.size() && !is_function && !is_no_data;
}

// ************************************************************************************
int CBuffer::gcd(int n, int m)
{
	while (n != m)
	{
		if (n > m)
			n -= m;
		else
			m -= n;
	}

	return n;
}

// EOF
