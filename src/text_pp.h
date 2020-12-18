#pragma once
// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2020-12-18
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

	void encode_word(vector<uint8_t>& v, uint32_t x);
	void encode_base(vector<uint8_t>& v, string str);
	void encode_plain(vector<uint8_t>& v, string str);
	void encode_number(vector<uint8_t>& v, string str);
	void encode_bars(vector<uint8_t>& v, uint32_t len);
	void encode_zero_run(vector<uint8_t>& v, uint32_t len);

	void decode_word(vector<uint8_t>& v_input, size_t& pos, vector<uint8_t>& v_output);
	void decode_number(vector<uint8_t>& v_input, size_t &pos, vector<uint8_t>& v_output);
	void decode_zero_run(vector<uint8_t>& v_input, size_t& pos, vector<uint8_t>& v_output);
	void decode_bars(vector<uint8_t>& v_input, size_t& pos, vector<uint8_t>& v_output);
	void decode_base(vector<uint8_t>& v, uint32_t c);

	token_t parse_word(vector<uint8_t>& v, uint32_t& pos, string& str);
	void parse_number(vector<uint8_t>& v, uint32_t& pos, string& str);
	void parse_bars(vector<uint8_t>& v, uint32_t& pos, string& str);
	void parse_zero_run(vector<uint8_t>& v, uint32_t& pos, string& str);

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
