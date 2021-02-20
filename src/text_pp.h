#pragma once
// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include <cstdio>
#include <algorithm>
#include <vector>
#include <iostream>
#include <string>
#include <tuple>
#include <list>
#include <unordered_map>
#include <map>
#include <cstdint>

using namespace std;

// ************************************************************************************
class CTextPreprocessing
{
	const uint32_t min_word_cnt = 16;
	const uint32_t min_word_len = 6;

	enum class token_t { nothing, word, number, bars, zero_run, base };

	unordered_map<string, uint32_t> m_dict;
	unordered_map<string, uint32_t> t_dict;
	
	bool word_symbol[256];

	vector<string> v_new_words;
	
	vector<string> v_dict;
	uint32_t dict_id;

	vector<pair<token_t, string>> v_tokens;

	// ************************************************
	template <size_t BASE, size_t RANGE>
	void encode_value(size_t x)
	{
		uint8_t s[32];
		int len = 0;

		for (; x; x /= RANGE)
		{
			size_t d = x % RANGE;
			s[len++] = (uint8_t)d;
		}
	}

	token_t get_token(vector<uint8_t> &v, uint32_t& pos, string& str);

	void encode_word(vector<uint8_t>& v, uint32_t x)
	{
		//	encode_value<128, 50>(x);

		if (x < 256)
		{
			v.emplace_back(255);
			v.emplace_back(x);

			return;
		}

		x -= 256;

		if (x < 256 * 256)
		{
			v.emplace_back(254);
			v.emplace_back(x >> 8);
			v.emplace_back(x & 0xff);

			return;
		}

		x -= 256 * 256;

		if (x < 256 * 256 * 256)
		{
			v.emplace_back(253);
			v.emplace_back(x >> 16);
			v.emplace_back((x >> 8) & 0xff);
			v.emplace_back(x & 0xff);

			return;
		}

		// !!! For future - consider even larger values
	}

	void encode_base(vector<uint8_t>& v, string &str)
	{
		if (str[0] == 'A')
			v.emplace_back(1);
		else if (str[0] == 'C')
			v.emplace_back(2);
		else if (str[0] == 'G')
			v.emplace_back(3);
		else if (str[0] == 'T')
			v.emplace_back(4);
	}

	void encode_plain(vector<uint8_t>& v, string &str)
	{
		for (auto c : str)
			v.emplace_back(c);
	}

	void encode_number(vector<uint8_t>& v, string &str)
	{
		if (str.size() > 15)
		{
			encode_plain(v, str);
			return;
		}

		// Range 128..227
		size_t x = 0;
		const size_t base = 100;

		for (auto c : str)
			x = x * 10 + (size_t)(c - '0');

		char t[16];
		int t_len = 0;

		for (; x; x /= base)
			t[t_len++] = (char)(x % base);

		for (int i = t_len - 1; i >= 0; --i)
			v.emplace_back(128 + t[i]);
	}

	void encode_bars(vector<uint8_t>& v, uint32_t len)
	{
		while (len)
		{
			int x = (len > 15) ? 15 : len;
			v.emplace_back(253 - x);				// 238..252
			len -= x;
		}
	}

	void encode_zero_run(vector<uint8_t>& v, uint32_t len)
	{
		while (len)
		{
			int x = (len > 10) ? 10 : len;
			v.emplace_back(238 - x);				// 228..237
			len -= x;
		}
	}

	void decode_word(vector<uint8_t>& v_input, size_t& pos, vector<uint8_t>& v_output)
	{
		int prefix = v_input[pos++];
		int code = 0;

		if (prefix == 255)
		{
			code = v_input[pos++];
		}
		else if (prefix == 254)
		{
			code = 256;

			code += (int)v_input[pos++] * 256;
			code += (int)v_input[pos++];
		}
		else if (prefix == 253)
		{
			code = 256 + 256 * 256;

			code += (int)v_input[pos++] * 256 * 256;
			code += (int)v_input[pos++] * 256;
			code += (int)v_input[pos++];
		}

		v_output.insert(v_output.end(), v_dict[code].begin(), v_dict[code].end());
	}

	void decode_number(vector<uint8_t>& v_input, size_t &pos, vector<uint8_t>& v_output)
	{
		// Range 128..227
		size_t x = 0;
		const size_t base = 100;

		while (pos < v_input.size())
		{
			int c = v_input[pos++];

			if (c < 128 || c > 227)
			{
				--pos;
				break;
			}

			x = x * base + (size_t)(c - 128);
		}

		string s = to_string(x);

		v_output.insert(v_output.end(), s.begin(), s.end());
	}

	void decode_zero_run(vector<uint8_t>& v_input, size_t& pos, vector<uint8_t>& v_output)
	{
		while (pos < v_input.size())
		{
			int c = v_input[pos++];

			if (c < 228 || c > 237)
			{
				--pos;
				break;
			}

			int len = 238 - c;

			for (int i = 0; i < len; ++i)
				v_output.emplace_back('0');
		}
	}

	void decode_bars(vector<uint8_t>& v_input, size_t& pos, vector<uint8_t>& v_output)
	{
		while (pos < v_input.size())
		{
			int c = v_input[pos++];

			if (c < 238 || c > 252)
			{
				--pos;
				break;
			}

			int len = 253 - c;

			for (int i = 0; i < len; ++i)
				v_output.emplace_back('|');
		}
	}

	void decode_base(vector<uint8_t>& v, uint32_t c)
	{
		v.emplace_back(" ACGT"[c]);
		v.emplace_back(':');
	}

	token_t parse_word(vector<uint8_t>& v, uint32_t& pos, string& str)
	{
		str.push_back(v[pos++]);

		if (pos == v.size())
			return token_t::word;

		if (v[pos] == ':' && (str[0] == 'A' || str[0] == 'C' || str[0] == 'G' || str[0] == 'T'))
		{
			str.push_back(v[pos]);
			++pos;

			return token_t::base;
		}

		while (pos < v.size())
		{
			auto c = v[pos++];

			if(word_symbol[c])
				str.push_back((char)c);
			else
			{
				--pos;
				break;
			}
		}

		return token_t::word;
	}

	void parse_number(vector<uint8_t>& v, uint32_t& pos, string& str)
	{
		for (; pos < v.size() && v[pos] >= '0' && v[pos] <= '9'; ++pos)
			str.push_back(v[pos]);
	}

	void parse_bars(vector<uint8_t>& v, uint32_t& pos, string& str)
	{
		for (; pos < v.size() && v[pos] == '|'; ++pos)
			str.push_back('|');
	}

	void parse_zero_run(vector<uint8_t>& v, uint32_t& pos, string& str)
	{
		for (; pos < v.size() && v[pos] == '0'; ++pos)
			str.push_back('0');
	}

	void update_dict(vector<uint8_t> &v_input);
	void store_dict(vector<uint8_t> &v_output);
	void load_dict(vector<uint8_t> &v_input, size_t &pos);
	void compress_part(vector<uint8_t> &v_input, vector<uint8_t>& v_output);
	void decompress_part(vector<uint8_t> &v_input, size_t &pos, vector<uint8_t>& v_output);

public:
	CTextPreprocessing();
	~CTextPreprocessing();

	void EncodeText(vector<uint8_t>& v_input, vector<uint8_t>& v_output);
	void DecodeText(vector<uint8_t>& v_input, vector<uint8_t>& v_output);
};

// EOF
