// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include "text_pp.h"

// ************************************************************************************
CTextPreprocessing::CTextPreprocessing()
{
	dict_id = 0;

	fill_n(word_symbol, 256, false);

	for(int c = 'A'; c <= 'Z'; ++c)
		word_symbol[c] = true;
	
	for(int c = 'a'; c <= 'z'; ++c)
		word_symbol[c] = true;

	for(int c = '0'; c <= '9'; ++c)
		word_symbol[c] = true;
	
	word_symbol['_'] = true;
	word_symbol['('] = true;
	word_symbol[')'] = true;
	word_symbol['&'] = true;
	word_symbol['/'] = true;
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

// EOF
