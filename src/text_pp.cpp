// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2020-12-18
// *******************************************************************************************

#include "text_pp.h"

// ************************************************************************************
CTextPreprocessing::CTextPreprocessing()
{
	dict_id = 0;
}

// ************************************************************************************
CTextPreprocessing::~CTextPreprocessing()
{
}

// ************************************************************************************
void CTextPreprocessing::EncodeText(vector<uint8_t>& v_input, vector<uint8_t>& v_output)
{
	v_output.clear();

	update_dict(v_input);
	store_dict(v_output);
	compress_part(v_input, v_output);
}

// ************************************************************************************
void CTextPreprocessing::DecodeText(vector<uint8_t>& v_input, vector<uint8_t>& v_output)
{
	size_t pos = 0;

	v_output.clear();
	load_dict(v_input, pos);

	decompress_part(v_input, pos, v_output);
}

// ************************************************
void CTextPreprocessing::update_dict(vector<uint8_t> &v_input)
{
	string str;
	token_t token;
	uint32_t pos = 0;

	v_new_words.clear();
	v_tokens.clear();

	while (pos < v_input.size())
	{
		token = get_token(v_input, pos, str);

		v_tokens.emplace_back(token, str);

		if (token == token_t::word)
		{
			auto p = m_dict.find(str);

			if (p == m_dict.end())
			{
				auto q = t_dict.find(str);
				if (q == t_dict.end())
					t_dict.emplace(str, 1);
				else
				{
					q->second++;
					if (q->second == min_word_cnt)
					{
						m_dict.emplace(q->first, dict_id++);
						v_new_words.emplace_back(q->first);
						t_dict.erase(q);
					}
				}
			}
		}
	}
}

// ************************************************
void CTextPreprocessing::store_dict(vector<uint8_t>& v_output)
{
	for (auto& s : v_new_words)
	{
		v_output.insert(v_output.end(), s.begin(), s.end());
		v_output.emplace_back('\n');
	}

	v_output.emplace_back(0);
}

// ************************************************
void CTextPreprocessing::load_dict(vector<uint8_t>& v_input, size_t& pos)
{
	string s;

	while (true)
	{
		uint8_t c = v_input[pos++];

		if (c == 0)
		{
			if (!s.empty())
				v_dict.emplace_back(s);
			s.clear();
			break;
		}
		else if (c == '\n')
		{
			if (!s.empty())
				v_dict.emplace_back(s);
			s.clear();
		}
		else
			s.push_back(c);
	}
}

// ************************************************
void CTextPreprocessing::compress_part(vector<uint8_t>& v_input, vector<uint8_t>& v_output)
{
	for (auto& t : v_tokens)
	{
		if (t.first == token_t::word)
		{
			auto p = m_dict.find(t.second);

			if (p != m_dict.end())
				encode_word(v_output, p->second);
			else
				encode_plain(v_output, t.second);
		}
		else if (t.first == token_t::base)
			encode_base(v_output, t.second);
		else if (t.first == token_t::number)
			encode_number(v_output, t.second);
		else if (t.first == token_t::bars)
			encode_bars(v_output, (uint32_t) t.second.size());
		else if (t.first == token_t::zero_run)
			encode_zero_run(v_output, (uint32_t) t.second.size());
		else
			encode_plain(v_output, t.second);
	}
}

// ************************************************
void CTextPreprocessing::decompress_part(vector<uint8_t>& v_input, size_t& pos, vector<uint8_t>& v_output)
{
	while (pos < v_input.size())
	{
		int c = v_input[pos++];

		if (c > 0 && c < 5)
			decode_base(v_output, c);
		else if (c < 128)
			v_output.emplace_back(c);
		else if (c < 228)
		{
			--pos;
			decode_number(v_input, pos, v_output);
		}
		else if (c < 238)
		{
			--pos;
			decode_zero_run(v_input, pos, v_output);
		}
		else if (c < 253)
		{
			--pos;
			decode_bars(v_input, pos, v_output);
		}
		else
		{
			--pos;
			decode_word(v_input, pos, v_output);
		}
	}
}

// ************************************************
CTextPreprocessing::token_t CTextPreprocessing::get_token(vector<uint8_t>& v, uint32_t& pos, string& str)
{
	auto c = v[pos];

	str.clear();

	if (c >= '1' && c <= '9')
	{
		parse_number(v, pos, str);

		return token_t::number;
	}
	else if (c == '0')
	{
		parse_zero_run(v, pos, str);

		return token_t::zero_run;
	}
	else if (c == '|')
	{
		parse_bars(v, pos, str);

		return token_t::bars;
	}
	else if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '_')
	{
		auto t = parse_word(v, pos, str);

		if (t == token_t::base)
			return t;

		if (str.size() >= min_word_len)
			return token_t::word;
		else
			return token_t::nothing;
	}

	str.push_back((char)c);
	++pos;

	return token_t::nothing;
}

// ************************************************
CTextPreprocessing::token_t CTextPreprocessing::parse_word(vector<uint8_t>& v, uint32_t& pos, string& str)
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

		if (c >= 'A' && c <= 'Z')
			str.push_back((char)c);
		else if (c >= 'a' && c <= 'z')
			str.push_back((char)c);
		else if (c >= '0' && c <= '9')
			str.push_back((char)c);
		else if (c == '_' || c == '(' || c == ')' || c == '&' || c == '/')
			str.push_back((char)c);
		else
		{
			--pos;
			break;
		}
	}

	return token_t::word;
}

// ************************************************
void CTextPreprocessing::parse_number(vector<uint8_t>& v, uint32_t& pos, string& str)
{
	for (; pos < v.size() && v[pos] >= '0' && v[pos] <= '9'; ++pos)
		str.push_back(v[pos]);
}

// ************************************************
void CTextPreprocessing::parse_bars(vector<uint8_t>& v, uint32_t& pos, string& str)
{
	for (; pos < v.size() && v[pos] == '|'; ++pos)
		str.push_back('|');
}

// ************************************************
void CTextPreprocessing::parse_zero_run(vector<uint8_t>& v, uint32_t& pos, string& str)
{
	for (; pos < v.size() && v[pos] == '0'; ++pos)
		str.push_back('0');
}

// ************************************************
void CTextPreprocessing::encode_word(vector<uint8_t>& v, uint32_t x)
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

// ************************************************
void CTextPreprocessing::encode_base(vector<uint8_t>& v, string str)
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

// ************************************************
void CTextPreprocessing::encode_plain(vector<uint8_t>& v, string str)
{
	for (auto c : str)
		v.emplace_back(c);
}

// ************************************************
void CTextPreprocessing::encode_number(vector<uint8_t>& v, string str)
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
		t[t_len++] = (char) (x % base);

	for (int i = t_len - 1; i >= 0; --i)
		v.emplace_back(128 + t[i]);
}

// ************************************************
void CTextPreprocessing::encode_bars(vector<uint8_t>& v, uint32_t len)
{
	while (len)
	{
		int x = (len > 15) ? 15 : len;
		v.emplace_back(253 - x);				// 238..252
		len -= x;
	}
}

// ************************************************
void CTextPreprocessing::encode_zero_run(vector<uint8_t>& v, uint32_t len)
{
	while (len)
	{
		int x = (len > 10) ? 10 : len;
		v.emplace_back(238 - x);				// 228..237
		len -= x;
	}
}

// ************************************************
void CTextPreprocessing::decode_word(vector<uint8_t>& v_input, size_t& pos, vector<uint8_t>& v_output)
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
		code = 256 + 256*256;

		code += (int)v_input[pos++] * 256 * 256;
		code += (int)v_input[pos++] * 256;
		code += (int)v_input[pos++];
	}

	v_output.insert(v_output.end(), v_dict[code].begin(), v_dict[code].end());
}

// ************************************************
void CTextPreprocessing::decode_number(vector<uint8_t>& v_input, size_t &pos, vector<uint8_t>& v_output)
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

// ************************************************
void CTextPreprocessing::decode_zero_run(vector<uint8_t>& v_input, size_t& pos, vector<uint8_t>& v_output)
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

// ************************************************
void CTextPreprocessing::decode_bars(vector<uint8_t>& v_input, size_t& pos, vector<uint8_t>& v_output)
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

// ************************************************
void CTextPreprocessing::decode_base(vector<uint8_t>& v, uint32_t c)
{
	v.emplace_back(" ACGT"[c]);
	v.emplace_back(':');
}

// EOF
