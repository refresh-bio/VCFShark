#pragma once
// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include <algorithm>
#include <vector>
#include <cstdint>
#include <cstring>
#include "defs.h"

using namespace std;

#define IGNORE_OFFSET_IN_FIRST_BLOCK

// ************************************************************************************
class CBuffer {
public:
	enum class buffer_t {none, flag, integer, real, text};

private:
	uint32_t max_size;
#ifndef IGNORE_OFFSET_IN_FIRST_BLOCK
	uint32_t offset;
#endif
	buffer_t type;

	vector<uint32_t> v_size;
	vector<uint8_t> v_data;

	bool is_function;
	bool is_no_data;		// buffer was set to empty vector

	function_data_item_t fun;

	uint32_t v_size_pos;
	uint32_t v_data_pos;

	uint32_t encode_var_int(char* p)
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

	uint32_t decode_var_int(uint8_t* p, uint32_t &val)
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
			val = (uint32_t)(((int)code) - 250);
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

			val = (uint32_t)-((int)val);

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

			val = (uint32_t)-((int)val);

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

			val = (uint32_t)-((int)val);

			return 5;
		}

		return 0;	// !!! Never should be here
	}

	uint32_t encode_float(char* p)
	{
		v_data.emplace_back(p[0]);
		v_data.emplace_back(p[1]);
		v_data.emplace_back(p[2]);
		v_data.emplace_back(p[3]);

		return 4;
	}

	uint32_t decode_float(uint8_t* p, float &val)
	{
		uint32_t x = 0;

		x += p[3];		x <<= 8;
		x += p[2];		x <<= 8;
		x += p[1];		x <<= 8;
		x += p[0];

		memcpy(&val, &x, 4);

		return 4;
	}

	void permute_integer_series_forward();
	void permute_integer_series_backward();
	void permute_float_series_forward();
	void permute_float_series_backward();

	int gcd(int n, int m)
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

public:
	CBuffer();
	~CBuffer();

	void SetMaxSize(uint32_t _max_size, uint32_t _offset);

	// Output buffer methods
	void WriteFlag(uint8_t f)
	{
		v_size.emplace_back(f);

		type = buffer_t::flag;
	}

	void WriteInt(char* p, uint32_t size)
	{
		v_size.emplace_back(size);

		v_data.insert(v_data.end(), p, p + 4 * size);

		type = buffer_t::integer;
	}

	void WriteIntVarSize(char* p, uint32_t size)
	{
		v_size.emplace_back(size);

		for (size_t i = 0; i < 4u * size; i += 4u)
			encode_var_int(p + i);
	}

	void WriteInt64(int64_t x)
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

	void WriteReal(char* p, uint32_t size)
	{
		v_size.emplace_back(size);

		v_data.insert(v_data.end(), p, p + 4 * size);

		type = buffer_t::real;
	}

	void WriteText(char* p, uint32_t size)
	{
		v_size.push_back(size);

		if (size)
			v_data.insert(v_data.end(), p, p + size);

		type = buffer_t::text;
	}

	void GetBuffer(vector<uint32_t>& _v_size, vector<uint8_t>& _v_data);
	
	bool IsFull(void)
	{
#ifndef IGNORE_OFFSET_IN_FIRST_BLOCK
		return v_data.size() + 4 * v_size.size() >= max_size + offset;
#else
		return v_data.size() + 4 * v_size.size() >= max_size;
#endif
	}

	// Input buffer methods
	void ReadFlag(uint8_t &flag)
	{
		if (v_size.empty())
		{
			flag = 0;

			return;
		}

		flag = (uint8_t)v_size[v_size_pos++];
	}

	void ReadInt(char* &p, uint32_t& size)
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

	void ReadIntVarSize(char* &p, uint32_t& size)
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

	void ReadInt64(int64_t &x)
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

	void ReadReal(char* &p, uint32_t& size)
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

	void ReadText(char* &p, uint32_t& size)
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

	void SetBuffer(vector<uint32_t>& _v_size, vector<uint8_t>& _v_data);

	void SetFunction(function_data_item_t& _fun);
	void FuncInt(char*& p, uint32_t& size, char* src_p, uint32_t src_size);
	void FuncReal(char*& p, uint32_t& size, char* src_p, uint32_t src_size);

	bool IsEmpty()
	{
		return v_size_pos >= v_size.size() && !is_function && !is_no_data;
	}
};

// EOF
