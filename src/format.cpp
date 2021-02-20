// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include "format.h"
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <limits>

//#define LOG_INFO

//#define USE_COMPRESSION_LEVEL

// *****************************************************************************************
CFormatCompress::CFormatCompress(string _desc, uint32_t _compression_level) : dict(ht_empty_key, 16, 0.6)		// assumption: 0x7fffffffu is unusual and should not appear
{
	desc = _desc;
	compression_level = _compression_level;

	no_samples = 1;
		
	vios_i = new CVectorIOStream(v_vios_i);
	vios_o = new CVectorIOStream(v_vios_o);

	rcd = new CRangeDecoder<CVectorIOStream>(*vios_i);
	rce = new CRangeEncoder<CVectorIOStream>(*vios_o);
}

// *****************************************************************************************
CFormatCompress::~CFormatCompress()
{
#ifdef LOG_INFO
	cerr << "***** " << desc << endl;
	cerr << "ctx_map_same        : " << to_string(ctx_map_same.get_size()) << endl;
	cerr << "ctx_map_known       : " << to_string(ctx_map_known.get_size()) << endl;
	cerr << "ctx_map_known2      : " << to_string(ctx_map_known2.get_size()) << endl;
	cerr << "ctx_map_plain       : " << to_string(ctx_map_plain.get_size()) << endl;
	cerr << "ctx_map_code        : " << to_string(ctx_map_code.get_size()) << endl;
	cerr << "ctx_map_entropy_type: " << to_string(ctx_map_entropy_type.get_size()) << endl;
#endif

	delete rcd;
	delete rce;

	delete vios_i;
	delete vios_o;
}

// *****************************************************************************************
void CFormatCompress::SetNoSamples(uint32_t _no_samples)
{
	no_samples = _no_samples;
}

// *****************************************************************************************
pair<CFormatCompress::info_t, uint32_t> CFormatCompress::determine_info_type(vector<uint32_t>& v_size)
{
	unordered_set<uint32_t> values;
	
	for (auto x : v_size)
	{
		values.insert(x);

		if (values.size() > 2)
			return make_pair(info_t::any, 0);
	}

	if (values.size() == 1)
	{
		if (*values.begin() == 0)
			return make_pair(info_t::zero, 0);
		if (*values.begin() == 1)
			return make_pair(info_t::one, 1);
		return make_pair(info_t::constant, *values.begin());
	}

	if (values.count(0) && values.count(1))
		return make_pair(info_t::zero_one, 1);
	
	if (values.count(0))
	{
		values.erase(0);
		return make_pair(info_t::zero_constant, *values.begin());
	}

	return make_pair(info_t::any, 0);
}

// *****************************************************************************************
template <unsigned SIZE> double CFormatCompress::entropy(vector<array<uint32_t, SIZE>> &vec)
{
	if (vec.size() == 1)
		return 0.0;

	double r = 0;

	uint32_t stats_size = 0;

	sort(vec.begin(), vec.end());

	vector<pair<uint32_t, uint32_t>> v_counts;

	auto p = vec.begin();
	v_counts.emplace_back(p->back(), 1);

	for (auto q = p + 1; q != vec.end(); ++q, ++p)
	{
		if (equal(p->begin(), p->begin() + SIZE - 1, q->begin()))
		{
			if (q->back() == v_counts.back().first)
				v_counts.back().second++;
			else
				v_counts.emplace_back(q->back(), 1);
		}
		else
		{
			++stats_size;
			r += entropy_0(v_counts);
			v_counts.clear();
			v_counts.emplace_back(q->back(), 1);
		}
	}

	++stats_size;
	r += entropy_0(v_counts);

	r += (double) stats_size * log2(stats_size);
	
	return r;
}

// *****************************************************************************************
void CFormatCompress::encode_ctx_type(int ctx_type)
{
	auto p_enc = find_rce_coder(ctx_map_entropy_type, 0);
	p_enc->Encode(ctx_type);
}

// *****************************************************************************************
int CFormatCompress::decode_ctx_type()
{
	auto p_dec = find_rcd_coder(ctx_map_entropy_type, 0);

	return p_dec->Decode();
}

// *****************************************************************************************
void CFormatCompress::encode_info_one(vector<uint8_t>& v_data, vector<uint8_t>& v_compressed)
{
	rce->Start();

	uint32_t* p_data = (uint32_t*)v_data.data();

	if (ctx_mode == 0)
	{
		array<double, 3> ent;

		{
			vector<array<uint32_t, 1>> vec;
			vec.reserve(v_data.size() / 4);
			for (uint32_t i = 0; i < v_data.size() / 4; ++i)
				vec.emplace_back(array<uint32_t, 1>{p_data[i]});
			ent[0] = entropy<1u>(vec);
		}

		{
			vector<array<uint32_t, 2>> vec;
			vec.reserve(v_data.size() / 4);
			for (uint32_t i = 1; i < v_data.size() / 4; ++i)
				vec.emplace_back(array<uint32_t, 2>{p_data[i - 1], p_data[i]});
			ent[1] = entropy<2u>(vec);
		}

		{
			vector<array<uint32_t, 3>> vec;
			vec.reserve(v_data.size() / 4);
			for (uint32_t i = 2; i < v_data.size() / 4; ++i)
				vec.emplace_back(array<uint32_t, 3>{p_data[i - 2], p_data[i - 1], p_data[i]});
			ent[2] = entropy<3u>(vec);
		}

		ctx_mode = (uint32_t) (min_element(ent.begin(), ent.end()) - ent.begin()) + 1;

		encode_ctx_type(ctx_mode);
	}
	
	context_t ctx = ~0ull;
	context_t ctx_mask;

	if (ctx_mode == 1)
		ctx_mask = 0xfffffull;
	else if (ctx_mode == 2)
		ctx_mask = 0xffffffffffull;
	else
		ctx_mask = 0xfffffffffffffffull;

	// Both floats and integers are treated as uint32_t as what we need is just to distinguish between different values
	uint32_t* q = p_data;
	auto p_enc = find_rce_coder(ctx_map_known, 0u);

	for (uint32_t i = 0; 4 * i < v_data.size(); ++i, ++q)
	{
		auto p = (*q == ht_empty_key) ? dict.local_end() : dict.find(*q);

		uint32_t code;

		ctx <<= 20;
		ctx &= ctx_mask;

		if (p == dict.local_end())
		{
			// Encode value plain
			p_enc->Encode(0);

			encode_plain(*q);

			code = (uint32_t) dict.size();
			if(code < max_dict_size && *q != ht_empty_key)
				dict.insert_fast(make_pair(*q, code));

			update_code_enc(ctx, code);
		}
		else
		{
			code = p->second;
			// Encode value as a code
			p_enc->Encode(1);

			encode_code(ctx, code);
		}

		ctx += code & 0xfffffull;
	}

	rce->End();

	v_compressed = move(v_vios_o);
}

// *****************************************************************************************
void CFormatCompress::decode_info_one(vector<uint8_t>& v_compressed, vector<uint8_t>& v_data)
{
	rcd->Start();

	if (ctx_mode == 0)
		ctx_mode = decode_ctx_type();

	context_t ctx = ~0ull;
	context_t ctx_mask;

	if (ctx_mode == 1)
		ctx_mask = 0xfffffull;
	else if (ctx_mode == 2)
		ctx_mask = 0xffffffffffull;
	else
		ctx_mask = 0xfffffffffffffffull;

	uint32_t max_i = (uint32_t) v_data.size();
	v_data.clear();

	auto p_dec = find_rcd_coder(ctx_map_known, 0u);

	// Both floats and integers are treated as uint32_t as what we need is just to distinguish between different values
	for (uint32_t i = 0; 4 * i < max_i; ++i)
	{
		uint32_t code;
		uint32_t val;

		ctx <<= 20;
		ctx &= ctx_mask;

		if (p_dec->Decode() == 0)
		{
			val = decode_plain();

			code = (uint32_t) dict_dec.size();
			if (code < max_dict_size && val != ht_empty_key)
				dict_dec.emplace_back(val);

			update_code_dec(ctx, code);
		}
		else
		{
			code = decode_code(ctx);
			val = dict_dec[code];
		}

		append_uint32(v_data, val);

		ctx += code & 0xfffffull;
	}

	rcd->End();
}

// *****************************************************************************************
void CFormatCompress::encode_info_constant(uint32_t val, vector<uint8_t>& v_data, vector<uint8_t>& v_compressed)
{
	rce->Start();

	uint32_t* p_data = (uint32_t*)v_data.data();

	uint32_t m = numeric_limits<uint32_t>::max();
	uint32_t s = val;
	uint32_t no_rows = (uint32_t) (v_data.size() / (4 * s));
	uint32_t no_items = (uint32_t) v_data.size() / 4;

	if (ctx_mode == 0)
	{
		array<double, 9> ent;
		// Entropy calculation for:
		// . . e - values
		// . d c - values
		// b a x - values
		//     p - pos

		// 0 - x
		{
			vector<array<uint32_t, 1>> vec;
			vec.reserve(no_items);
			for (uint32_t i = 0; i < no_items; ++i)
				vec.emplace_back(array<uint32_t, 1>{p_data[i]});
			ent[0] = entropy<1u>(vec);
		}

		// 1 - a, x
		{
			vector<array<uint32_t, 2>> vec;
			vec.reserve(no_rows * s);

			for (uint32_t i = 0; i < no_rows; ++i)
				for (uint32_t j = 0; j < s; ++j)
					vec.emplace_back(array<uint32_t, 2>{j == 0 ? m : p_data[i * s + j - 1], p_data[i * s + j]});
			ent[1] = entropy<2u>(vec);
		}

		// 2 - b, a, x
		{
			vector<array<uint32_t, 3>> vec;
			vec.reserve(no_rows * s);

			for (uint32_t i = 0; i < no_rows; ++i)
				for (uint32_t j = 0; j < s; ++j)
					vec.emplace_back(array<uint32_t, 3>{j < 2 ? m : p_data[i * s + j - 2], j == 0 ? m : p_data[i * s + j - 1], p_data[i * s + j]});
			ent[2] = entropy<3u>(vec);
		}

		// 3 - c, x
		{
			vector<array<uint32_t, 2>> vec;
			vec.reserve(no_rows * s);

			for (uint32_t i = 0; i < no_rows; ++i)
				for (uint32_t j = 0; j < s; ++j)
					vec.emplace_back(array<uint32_t, 2>{i == 0 ? m : p_data[i * s + j - s], p_data[i * s + j]});
			ent[3] = entropy<2u>(vec);
		}

		// 4 - e, c, x
		{
			vector<array<uint32_t, 3>> vec;
			vec.reserve(no_rows * s);

			for (uint32_t i = 0; i < no_rows; ++i)
				for (uint32_t j = 0; j < s; ++j)
					vec.emplace_back(array<uint32_t, 3>{i < 2 ? m : p_data[i * s + j - 2 * s], i == 0 ? m : p_data[i * s + j - s], p_data[i * s + j]});
			ent[4] = entropy<3u>(vec);
		}

		// 5 - a, c, x
		{
			vector<array<uint32_t, 3>> vec;
			vec.reserve(no_rows * s);

			for (uint32_t i = 0; i < no_rows; ++i)
				for (uint32_t j = 0; j < s; ++j)
					vec.emplace_back(array<uint32_t, 3>{j == 0 ? m : p_data[i * s + j - 1], i == 0 ? m : p_data[i * s + j - s], p_data[i * s + j]});
			ent[5] = entropy<3u>(vec);
		}

		// 6 - p, a, x
		{
			vector<array<uint32_t, 3>> vec;
			vec.reserve(no_rows * s);

			for (uint32_t i = 0; i < no_rows; ++i)
				for (uint32_t j = 0; j < s; ++j)
					vec.emplace_back(array<uint32_t, 3>{j, j == 0 ? m : p_data[i * s + j - 1], p_data[i * s + j]});
			ent[6] = entropy<3u>(vec);
		}

		// 7 - p, c, x
		{
			vector<array<uint32_t, 3>> vec;
			vec.reserve(no_rows * s);

			for (uint32_t i = 0; i < no_rows; ++i)
				for (uint32_t j = 0; j < s; ++j)
					vec.emplace_back(array<uint32_t, 3>{j, i == 0 ? m : p_data[i * s + j - s], p_data[i * s + j]});
			ent[7] = entropy<3u>(vec);
		}

		// 8 - p, x
		{
			vector<array<uint32_t, 2>> vec;
			vec.reserve(no_rows * s);

			for (uint32_t i = 0; i < no_rows; ++i)
				for (uint32_t j = 0; j < s; ++j)
					vec.emplace_back(array<uint32_t, 2>{j, p_data[i * s + j]});
			ent[8] = entropy<2u>(vec);
		}

		ctx_mode = (uint32_t) (min_element(ent.begin(), ent.end()) - ent.begin()) + 1;

		encode_ctx_type(ctx_mode);
	}

	context_t ctx = 0;
	const context_t ctx_small_mask = 0xfffffull;

	// Both floats and integers are treated as uint32_t as what we need is just to distinguish between different values
	uint32_t* q = p_data;

	vector<uint32_t> v_codes(no_items);

	auto p_enc = find_rce_coder(ctx_map_known, 0u);

	for (uint32_t i = 0; i < no_rows; ++i)
		for (uint32_t j = 0; j < s; ++j, ++q)
		{
			auto p = (*q == ht_empty_key) ? dict.local_end() : dict.find(*q);

			uint32_t code;

			switch (ctx_mode)
			{
			case 1: 
				ctx = 0;
				break;
			case 2: 
				ctx = ((j == 0) ? m : v_codes[i * s + j - 1]) & ctx_small_mask;
				break;
			case 3: 
				ctx = (((j < 2) ? m : v_codes[i * s + j - 2]) & ctx_small_mask) << 20;
				ctx += ((j == 0) ? m : v_codes[i * s + j - 1]) & ctx_small_mask;
				break;
			case 4:
				ctx = ((i == 0) ? m : v_codes[i * s + j - s]) & ctx_small_mask;
				break;
			case 5:
				ctx = (((i < 2) ? m : v_codes[i * s + j - 2 * s]) & ctx_small_mask) << 20;
				ctx += ((i == 0) ? m : v_codes[i * s + j - s]) & ctx_small_mask;
				break;
			case 6:
				ctx = (((j == 0) ? m : v_codes[i * s + j - 1]) & ctx_small_mask) << 20;
				ctx += ((i == 0) ? m : v_codes[i * s + j - s]) & ctx_small_mask;
				break;
			case 7:
				ctx = ((j == 0) ? m : v_codes[i * s + j - 1]) & ctx_small_mask;
				ctx += (j & ctx_small_mask) << 20;
				break;
			case 8:
				ctx = ((i == 0) ? m : v_codes[i * s + j - s]) & ctx_small_mask;
				ctx += (j & ctx_small_mask) << 20;
				break;
			case 9:
				ctx += j & ctx_small_mask;
				break;
			}

			if (p == dict.local_end())
			{
				// Encode value plain
				p_enc->Encode(0);

				encode_plain(*q);

				code = (uint32_t) dict.size();
				if(code < max_dict_size && *q != ht_empty_key)
					dict.insert_fast(make_pair(*q, code));

				update_code_enc(ctx, code);
			}
			else
			{
				code = p->second;
				// Encode value as a code
				p_enc->Encode(1);

				encode_code(ctx, code);
			}

			v_codes[i * s + j] = code;
		}

	rce->End();

	v_compressed = move(v_vios_o);
}

// *****************************************************************************************
void CFormatCompress::decode_info_constant(uint32_t val, vector<uint8_t>& v_compressed, vector<uint8_t>& v_data)
{
	rcd->Start();

	uint32_t m = numeric_limits<uint32_t>::max();
	uint32_t s = val;
	uint32_t no_rows = (uint32_t) (v_data.size() / (4 * s));
	uint32_t no_items = (uint32_t) v_data.size() / 4;

	if (ctx_mode == 0)
		ctx_mode = decode_ctx_type();

	context_t ctx = 0;
	context_t ctx_small_mask = 0xfffffull;

	v_data.clear();

	// Both floats and integers are treated as uint32_t as what we need is just to distinguish between different values
	vector<uint32_t> v_codes(no_items);

	auto p_dec = find_rcd_coder(ctx_map_known, 0u);

	for (uint32_t i = 0; i < no_rows; ++i)
		for (uint32_t j = 0; j < s; ++j)
		{
			uint32_t code;

			switch (ctx_mode)
			{
			case 1:
				ctx = 0;
				break;
			case 2:
				ctx = ((j == 0) ? m : v_codes[i * s + j - 1]) & ctx_small_mask;
				break;
			case 3:
				ctx = (((j < 2) ? m : v_codes[i * s + j - 2]) & ctx_small_mask) << 20;
				ctx += ((j == 0) ? m : v_codes[i * s + j - 1]) & ctx_small_mask;
				break;
			case 4:
				ctx = ((i == 0) ? m : v_codes[i * s + j - s]) & ctx_small_mask;
				break;
			case 5:
				ctx = (((i < 2) ? m : v_codes[i * s + j - 2 * s]) & ctx_small_mask) << 20;
				ctx += ((i == 0) ? m : v_codes[i * s + j - s]) & ctx_small_mask;
				break;
			case 6:
				ctx = (((j == 0) ? m : v_codes[i * s + j - 1]) & ctx_small_mask) << 20;
				ctx += ((i == 0) ? m : v_codes[i * s + j - s]) & ctx_small_mask;
				break;
			case 7:
				ctx = ((j == 0) ? m : v_codes[i * s + j - 1]) & ctx_small_mask;
				ctx += (j & ctx_small_mask) << 20;
				break;
			case 8:
				ctx = ((i == 0) ? m : v_codes[i * s + j - s]) & ctx_small_mask;
				ctx += (j & ctx_small_mask) << 20;
				break;
			case 9:
				ctx += j & ctx_small_mask;
				break;
			}

			if (p_dec->Decode() == 0)
			{
				val = decode_plain();

				code = (uint32_t) dict_dec.size();
				if (code < max_dict_size && val != ht_empty_key)
					dict_dec.emplace_back(val);

				update_code_dec(ctx, code);
			}
			else
			{
				code = decode_code(ctx);
				val = dict_dec[code];
			}

			append_uint32(v_data, val);

			v_codes[i * s + j] = code;
		}

	rcd->End();
}

// *****************************************************************************************
void CFormatCompress::encode_info_any(vector<uint32_t>& v_size, vector<uint8_t>& v_data, vector<uint8_t>& v_compressed)
{
	encode_info_one(v_data, v_compressed);
}

// *****************************************************************************************
void CFormatCompress::decode_info_any(vector<uint32_t>& v_size, vector<uint8_t>& v_compressed, vector<uint8_t>& v_data)
{
	decode_info_one(v_compressed, v_data);
}

// *****************************************************************************************
void CFormatCompress::encode_plain(uint32_t x)
{
	for (int k = 0; k < 4; ++k)
	{
		auto p_enc = find_rce_coder(ctx_map_plain, k);
		p_enc->Encode((x >> (8 * k)) & 0xff);
	}
}

// *****************************************************************************************
uint32_t CFormatCompress::decode_plain()
{
	uint32_t x = 0;

	for (int k = 0; k < 4; ++k)
	{
		auto p_dec = find_rcd_coder(ctx_map_plain, k);
		x += ((uint32_t) p_dec->Decode()) << (8 * k);
	}

	return x;
}

// *****************************************************************************************
void CFormatCompress::encode_code(context_t ctx, uint32_t code)
{
	uint32_t dict_size = (uint32_t) dict.size();

	if (dict_size > 256 * 256 * 256)
	{
		auto p_enc = find_rce_coder(ctx_map_code, ctx + (3ull << 62));
		p_enc->Encode(code >> 24);
	}

	if (dict_size > 256 * 256)
	{
		auto p_enc = find_rce_coder(ctx_map_code, ctx + (2ull << 62));
		p_enc->Encode((code >> 16) & 0xff);
#ifdef USE_COMPRESSION_LEVEL
		if(compression_level == 3)
#endif
			ctx += code & 0xf0000ull;
	}

	if (dict_size > 256)
	{
		auto p_enc = find_rce_coder(ctx_map_code, ctx + (1ull << 62));
		p_enc->Encode((code >> 8) & 0xff);

#ifdef USE_COMPRESSION_LEVEL
		if (compression_level >= 2)
			ctx += code & 0xff00ull;
		else
			ctx += code & 0x7f00ull;
#else
		ctx += code & 0xff00ull;
#endif
	}

	auto p_enc = find_rce_coder(ctx_map_code, ctx + (0ull << 62));
	p_enc->Encode(code & 0xff);
}

// *****************************************************************************************
uint32_t CFormatCompress::decode_code(context_t ctx)
{
	uint32_t dict_size = (uint32_t) dict_dec.size();
	uint32_t code = 0;

	if (dict_size > 256 * 256 * 256)
	{
		auto p_dec = find_rcd_coder(ctx_map_code, ctx + (3ull << 62));
		code = ((uint32_t) p_dec->Decode()) << 24;
	}

	if (dict_size > 256 * 256)
	{
		auto p_dec = find_rcd_coder(ctx_map_code, ctx + (2ull << 62));
		code += ((uint32_t) p_dec->Decode()) << 16;
#ifdef USE_COMPRESSION_LEVEL
		if(compression_level == 3)
#endif
			ctx += code & 0xf0000ull;
	}

	if (dict_size > 256)
	{
		auto p_dec = find_rcd_coder(ctx_map_code, ctx + (1ull << 62));
		code += ((uint32_t) p_dec->Decode()) << 8;
#ifdef USE_COMPRESSION_LEVEL
		if(compression_level >= 2)
			ctx += code & 0xff00ull;
		else
			ctx += code & 0x7f00ull;
#else
		ctx += code & 0xff00ull;
#endif
	}

	auto p_dec = find_rcd_coder(ctx_map_code, ctx + (0ull << 62));
	code += (uint32_t) p_dec->Decode();

	return code;
}

// *****************************************************************************************
void CFormatCompress::update_code_enc(context_t ctx, uint32_t code)
{
	uint32_t dict_size = (uint32_t) dict.size();

	if (dict_size > 256 * 256 * 256)
	{
		auto p_enc = find_rce_coder(ctx_map_code, ctx + (3ull << 62));
		p_enc->Update(code >> 24);
	}

	if (dict_size > 256 * 256)
	{
		auto p_enc = find_rce_coder(ctx_map_code, ctx + (2ull << 62));
		p_enc->Update((code >> 16) & 0xff);
#ifdef USE_COMPRESSION_LEVEL
		if(compression_level == 3)
#endif
			ctx += code & 0xf0000ull;
	}

	if (dict_size > 256)
	{
		auto p_enc = find_rce_coder(ctx_map_code, ctx + (1ull << 62));
		p_enc->Update((code >> 8) & 0xff);
#ifdef USE_COMPRESSION_LEVEL
		if (compression_level >= 2)
			ctx += code & 0xff00ull;
		else
			ctx += code & 0x7f00ull;
#else
		ctx += code & 0xff00ull;
#endif
	}

	auto p_enc = find_rce_coder(ctx_map_code, ctx + (0ull << 62));
	p_enc->Update(code & 0xff);
}

// *****************************************************************************************
void CFormatCompress::update_code_dec(context_t ctx, uint32_t code)
{
	uint32_t dict_size = (uint32_t) dict_dec.size();

	if (dict_size > 256 * 256 * 256)
	{
		auto p_dec = find_rcd_coder(ctx_map_code, ctx + (3ull << 62));
		p_dec->Update(code >> 24);
	}

	if (dict_size > 256 * 256)
	{
		auto p_dec = find_rcd_coder(ctx_map_code, ctx + (2ull << 62));
		p_dec->Update((code >> 16) & 0xff);
#ifdef USE_COMPRESSION_LEVEL
		if (compression_level == 3)
#endif
			ctx += code & 0xf0000ull;
	}

	if (dict_size > 256)
	{
		auto p_dec = find_rcd_coder(ctx_map_code, ctx + (1ull << 62));
		p_dec->Update((code >> 8) & 0xff);
#ifdef USE_COMPRESSION_LEVEL
		if (compression_level >= 2)
			ctx += code & 0xff00ull;
		else
			ctx += code & 0x7f00ull;
#else
		ctx += code & 0xff00ull;
#endif
	}

	auto p_dec = find_rcd_coder(ctx_map_code, ctx + (0ull << 62));
	p_dec->Update(code & 0xff);
}

// *****************************************************************************************
void CFormatCompress::encode_format_one(vector<uint32_t>& v_size, vector<uint8_t>& v_data, vector<uint8_t>& v_compressed)
{
	rce->Start();

	// Both floats and integers are treated as uint32_t as what we need is just to distinguish between different values
	uint32_t* q = (uint32_t*)v_data.data();

	array<vector<uint64_t>, 3> v_codes = { vector<uint64_t>(no_samples, 0), vector<uint64_t>(no_samples, 0), vector<uint64_t>(no_samples, 0) };

	auto p_enc = find_rce_coder(ctx_map_known, 0u);
	ctx_map_256_16_1_t::value_type p_enc_plain[4] = { find_rce_coder(ctx_map_plain, 0), find_rce_coder(ctx_map_plain, 1), find_rce_coder(ctx_map_plain, 2), find_rce_coder(ctx_map_plain, 3) };

	for (uint32_t i = 0; i < v_size.size(); ++i)
	{
		for (uint32_t j = 0; j < no_samples; ++j, ++q)
		{
			auto p = (*q == ht_empty_key) ? dict.local_end() : dict.find(*q);

			uint32_t code;

			if (p == dict.local_end())
			{
				// Encode value plain
				p_enc->Encode(0);

				for (int k = 0; k < 4; ++k)
					p_enc_plain[k]->Encode((*q >> (8 * k)) & 0xff);

				code = (uint32_t) dict.size();
				if (code < max_dict_size && *q != ht_empty_key)
					dict.insert_fast(make_pair(*q, code));
			}
			else
			{
				// Encode value as a code
				p_enc->Encode(1);

				context_t ctx = 0;

				if (i > 0)
					ctx += (v_codes[(i-1) % 3][j] & 0xfffffull) << 24;
				else
					ctx += 0xfffffull << 24;
				if (j > 0)
					ctx += (v_codes[i % 3][j - 1] & 0xfffffull) << 44;
				else
					ctx += 0xfffffull << 44;

				uint32_t dict_size = (uint32_t) dict.size();
				code = p->second;

				if (dict_size > 256 * 256 * 256)
				{
					auto p_enc = find_rce_coder(ctx_map_code, ctx + (3ull << 20));
					p_enc->Encode(code >> 24);
				}
				if (dict_size > 256 * 256)
				{
					auto p_enc = find_rce_coder(ctx_map_code, ctx + (2ull << 20));
					p_enc->Encode((code >> 16) & 0xff);
#ifdef USE_COMPRESSION_LEVEL
					if (compression_level == 3)
#endif
						ctx += code & 0xf0000;
				}
				if (dict_size > 256)
				{
					auto p_enc = find_rce_coder(ctx_map_code, ctx + (1ull << 20));
					p_enc->Encode((code >> 8) & 0xff);
#ifdef USE_COMPRESSION_LEVEL
					if (compression_level >= 2)
						ctx += code & 0xff00;
					else
						ctx += code & 0x7f00;
#else
					ctx += code & 0xff00;
#endif
				}
				auto p_enc = find_rce_coder(ctx_map_code, ctx + (0ull << 20));
				p_enc->Encode(code & 0xff);
			}

			v_codes[i % 3][j] = code;
		}
	}

	rce->End();

	v_compressed = move(v_vios_o);
}

// *****************************************************************************************
void CFormatCompress::decode_format_one(vector<uint32_t>& v_size, vector<uint8_t>& v_data, vector<uint8_t>& v_compressed)
{
	rcd->Start();

	// Both floats and integers are treated as uint32_t as what we need is just to distinguish between different values
	v_data.clear();
//	uint32_t* q = (uint32_t*)v_data.data();

	array<vector<uint64_t>, 3> v_codes = { vector<uint64_t>(no_samples, 0), vector<uint64_t>(no_samples, 0), vector<uint64_t>(no_samples, 0) };

	auto p_dec = find_rcd_coder(ctx_map_known, 0u);
	ctx_map_256_16_1_t::value_type p_dec_plain[4] = { find_rcd_coder(ctx_map_plain, 0), find_rcd_coder(ctx_map_plain, 1), find_rcd_coder(ctx_map_plain, 2), find_rcd_coder(ctx_map_plain, 3) };

	for (uint32_t i = 0; i < v_size.size(); ++i)
	{
		for (uint32_t j = 0; j < no_samples; ++j) //, ++q)
		{
			uint32_t code;
			uint32_t val;

			if (p_dec->Decode() == 0)
			{
				val = 0;

				for (int k = 0; k < 4; ++k)
				{
					uint32_t x = p_dec_plain[k]->Decode();
					val += x << (8 * k);
				}

				code = (uint32_t) dict_dec.size();
				if (code < max_dict_size && val != ht_empty_key)
					dict_dec.emplace_back(val);
			}
			else
			{
				context_t ctx = 0;

				if (i > 0)
					ctx += (v_codes[(i - 1) % 3][j] & 0xfffffull) << 24;
				else
					ctx += 0xfffffull << 24;
				if (j > 0)
					ctx += (v_codes[i % 3][j - 1] & 0xfffffull) << 44;
				else
					ctx += 0xfffffull << 44;

				uint32_t dict_size = (uint32_t) dict_dec.size();
				code = 0;

				if (dict_size > 256 * 256 * 256)
				{
					auto p_dec = find_rcd_coder(ctx_map_code, ctx + (3ull << 20));
					code += ((uint32_t) p_dec->Decode()) << 24;
				}
				if (dict_size > 256 * 256)
				{
					auto p_enc = find_rcd_coder(ctx_map_code, ctx + (2ull << 20));
					code += ((uint32_t) p_enc->Decode()) << 16;
#ifdef USE_COMPRESSION_LEVEL
					if (compression_level == 3)
#endif
						ctx += code & 0xf0000;
				}
				if (dict_size > 256)
				{
					auto p_dec = find_rcd_coder(ctx_map_code, ctx + (1ull << 20));
					code += ((uint32_t) p_dec->Decode()) << 8;
#ifdef USE_COMPRESSION_LEVEL
					if (compression_level >= 2)
						ctx += code & 0xff00;
					else
						ctx += code & 0x7f00;
#else
					ctx += code & 0xff00;
#endif
				}
				auto p_dec = find_rcd_coder(ctx_map_code, ctx + (0ull << 20));
				code += (uint32_t) p_dec->Decode();

				val = dict_dec[code];
			}

			append_uint32(v_data, val);

			v_codes[i % 3][j] = code;
		}
	}

	rcd->End();
}

// *****************************************************************************************
void CFormatCompress::encode_format_many(vector<uint32_t>& v_size, vector<uint8_t>& v_data, vector<uint8_t>& v_compressed)
{
	rce->Start();

	// Both floats and integers are treated as uint32_t as what we need is just to distinguish between different values
	uint32_t* p_data = (uint32_t*)v_data.data();

	uint32_t* cur_line = p_data;

	int cur_items_per_sample = 0;
	int prev_items_per_sample = 0;

	auto p_enc = find_rce_coder(ctx_map_same, 0u);
	auto p_enc2 = find_rce_coder(ctx_map_known2, 0u);
	ctx_map_256_16_1_t::value_type p_enc_plain[4] = { find_rce_coder(ctx_map_plain, 0), find_rce_coder(ctx_map_plain, 1), find_rce_coder(ctx_map_plain, 2), find_rce_coder(ctx_map_plain, 3) };

	for (uint32_t i = 0; i < v_size.size(); ++i)
	{
		int c_size = v_size[i];
		cur_items_per_sample = c_size / max(1u, no_samples);
		int entry_bytes = cur_items_per_sample * 4;

		if(c_size)
			for (uint32_t j = 0; j < max(1u, no_samples); ++j)
			{
				if (i > 0 && prev_items_per_sample == cur_items_per_sample)
				{
					if (memcmp(cur_line, cur_line - prev_items_per_sample, entry_bytes) == 0)
					{
						p_enc->Encode(1);
						cur_line += cur_items_per_sample;
						continue;
					}

					p_enc->Encode(0);
				}
			
				context_t ctx = 0u;

				for (int k = 0; k < cur_items_per_sample; ++k)
				{
					uint32_t x = cur_line[k];
					auto p = (x == ht_empty_key) ? dict.local_end() : dict.find(x);

					ctx &= 0x1fffffffffffull;
					ctx += ((uint64_t)k) << 58;

					if (k == 0)
					{
						ctx += ((uint64_t)j & 0x1fffull) << 45;
					}	

					uint32_t code;

					if (p == dict.local_end())
					{
						// Encode value in plain
						p_enc2->Encode(0);

						for (int a = 0; a < 4; ++a)
							p_enc_plain[a]->Encode((x >> (8 * a)) & 0xff);

						code = (uint32_t) dict.size();
						if (code < max_dict_size && x != ht_empty_key)
							dict.insert_fast(make_pair(x, code));

						ctx += code;

						ctx <<= 15;
					}
					else
					{
						p_enc2->Encode(1);
						uint32_t dict_size = (uint32_t) dict.size();
						code = p->second;

						if (dict_size > 256 * 256 * 256)
						{
							auto p_enc = find_rce_coder(ctx_map_code, ctx + (4ull << 61));
							p_enc->Encode(code >> 24);
						}
						if (dict_size > 256 * 256)
						{
							auto p_enc = find_rce_coder(ctx_map_code, ctx + (3ull << 61));
							p_enc->Encode((code >> 16) & 0xff);
						}
						if (dict_size > 256)
						{
							auto p_enc = find_rce_coder(ctx_map_code, ctx + (2ull << 61));
							p_enc->Encode((code >> 8) & 0xff);
#ifdef USE_COMPRESSION_LEVEL
							if (compression_level >= 2)
								ctx += code & 0x7f00;
							else
								ctx += code & 0x1f00;
#else
							ctx += code & 0x7f00;
#endif
						}
						auto p_enc = find_rce_coder(ctx_map_code, ctx + (1ull << 61));
						p_enc->Encode(code & 0xff);

						ctx += code & 0xff;

						ctx <<= 15;
					}
				}

				cur_line += cur_items_per_sample;
			}

		prev_items_per_sample = cur_items_per_sample;
	} 

	rce->End();

	v_compressed = move(v_vios_o);
}

// *****************************************************************************************
void CFormatCompress::decode_format_many(vector<uint32_t>& v_size, vector<uint8_t>& v_data, vector<uint8_t>& v_compressed)
{
	rcd->Start();

	// Both floats and integers are treated as uint32_t as what we need is just to distinguish between different values
	v_data.clear();
	v_data.reserve(accumulate(v_size.begin(), v_size.end(), 0u) * 4);

	uint8_t* p_data = (uint8_t*)v_data.data();
	uint32_t* cur_line = (uint32_t*) p_data;

	int cur_items_per_sample = 0;
	int prev_items_per_sample = 0;

/*#ifdef LOG_INFO
	cout << "Start decode_format_many\n";
#endif*/

	auto p_dec = find_rcd_coder(ctx_map_same, 0u);
	auto p_dec2 = find_rcd_coder(ctx_map_known2, 0u);
	ctx_map_256_16_1_t::value_type p_dec_plain[4] = { find_rcd_coder(ctx_map_plain, 0), find_rcd_coder(ctx_map_plain, 1),find_rcd_coder(ctx_map_plain, 2), find_rcd_coder(ctx_map_plain, 3) };

	for (uint32_t i = 0; i < v_size.size(); ++i)
	{
		int c_size = v_size[i];
		cur_items_per_sample = c_size / max(1u, no_samples);
		int entry_bytes = cur_items_per_sample * 4;

		if(c_size)
			for (uint32_t j = 0; j < max(1u, no_samples); ++j)
			{
				if (i > 0 && prev_items_per_sample == cur_items_per_sample)
				{
					if (p_dec->Decode() == 1)
					{
						int eb4 = entry_bytes / 4;
						v_data.resize(v_data.size() + entry_bytes);
						for (int k = 0; k < entry_bytes / 4; ++k)
						{
							*cur_line = *(cur_line - eb4);
							++cur_line;
						}
						continue;
					}
				}

				context_t ctx = 0u;

				for (int k = 0; k < cur_items_per_sample; ++k)
				{
					uint32_t val = 0;
					uint32_t code = 0;

					ctx &= 0x1fffffffffffull;
					ctx += ((uint64_t)k) << 58;

					if (k == 0)
					{
						ctx += ((uint64_t)j & 0x1fffull) << 45;
					}

					if (p_dec2->Decode() == 0)
					{
						for (int a = 0; a < 4; ++a)
							val += ((uint32_t) p_dec_plain[a]->Decode()) << (8 * a);

						code = (uint32_t) dict_dec.size();
						if (code < max_dict_size && val != ht_empty_key)
							dict_dec.emplace_back(val);

						ctx += code;
						ctx <<= 15;
					}
					else
					{
						uint32_t dict_size = (uint32_t) dict_dec.size();
						code = 0;

						if (dict_size > 256 * 256 * 256)
						{
							auto p_dec = find_rcd_coder(ctx_map_code, ctx + (4ull << 61));
							code += ((uint32_t)p_dec->Decode()) << 24;
						}
						if (dict_size > 256 * 256)
						{
							auto p_enc = find_rcd_coder(ctx_map_code, ctx + (3ull << 61));
							code += ((uint32_t)p_enc->Decode()) << 16;
						}
						if (dict_size > 256)
						{
							auto p_dec = find_rcd_coder(ctx_map_code, ctx + (2ull << 61));
							code += ((uint32_t)p_dec->Decode()) << 8;
#ifdef USE_COMPRESSION_LEVEL
							if (compression_level >= 2)
								ctx += code & 0x7f00;
							else
								ctx += code & 0x1f00;
#else
							ctx += code & 0x7f00;
#endif
						}
						auto p_dec = find_rcd_coder(ctx_map_code, ctx + (1ull << 61));
						code += (uint32_t)p_dec->Decode();

						val = dict_dec[code];

						ctx += code & 0xff;
						ctx <<= 15;
					}

					append_uint32(v_data, val);
				}

				cur_line += cur_items_per_sample;
			}

		prev_items_per_sample = cur_items_per_sample;
	} 

	rcd->End();
}

// *****************************************************************************************
void CFormatCompress::EncodeFormat(vector<uint32_t> &v_size, vector<uint8_t>& v_data, vector<uint8_t>& v_compressed)
{
	if (v_data.empty())
	{
		v_compressed.clear();
		return;
	}

	bool one = true;

	for(auto x : v_size)
		if (x != no_samples)
		{
			one = false;
			break;
		}

	if (one)
		encode_format_one(v_size, v_data, v_compressed);
	else
		encode_format_many(v_size, v_data, v_compressed);
}

// *****************************************************************************************
void CFormatCompress::DecodeFormat(vector<uint32_t>& v_size, vector<uint8_t>& v_compressed, vector<uint8_t>& v_data)
{
	if (v_compressed.empty())
	{
		v_data.clear();
		return;
	}

	bool one = true;

	for (auto x : v_size)
		if (x != no_samples)
		{
			one = false;
			break;
		}

	v_vios_i = move(v_compressed);
	vios_i->RestartRead();

	if (one)
		decode_format_one(v_size, v_data, v_compressed);
	else
		decode_format_many(v_size, v_data, v_compressed);
}

// *****************************************************************************************
void CFormatCompress::EncodeInfo(vector<uint32_t> &v_size, vector<uint8_t>& v_data, vector<uint8_t>& v_compressed)
{
	if (v_data.empty())
	{
		v_compressed.clear();
		return;
	}

	auto act_type = determine_info_type(v_size);

	if (act_type.first != type.first)
	{
		type = act_type;
		ctx_mode = 0;
	}

	if (type.first == info_t::zero)
		return;				// Never should be here
	if (type.first == info_t::one)
		encode_info_one(v_data, v_compressed);
	else if (type.first == info_t::zero_one)
		encode_info_one(v_data, v_compressed);
	else if (type.first == info_t::constant)
		encode_info_constant(type.second, v_data, v_compressed);
	else if (type.first == info_t::zero_constant)
		encode_info_constant(type.second, v_data, v_compressed);
	else
		encode_info_any(v_size, v_data, v_compressed);
}

// *****************************************************************************************
void CFormatCompress::DecodeInfo(vector<uint32_t>& v_size, vector<uint8_t>& v_compressed, vector<uint8_t>& v_data)
{
	if (v_compressed.empty())
	{
		v_data.clear();
		return;
	}

	auto act_type = determine_info_type(v_size);

	if (act_type.first != type.first)
	{
		type = act_type;
		ctx_mode = 0;
	}

	v_vios_i = move(v_compressed);
	vios_i->RestartRead();

	if (type.first == info_t::zero)
		return;				// Never should be here
	if (type.first == info_t::one)
		decode_info_one(v_compressed, v_data);
	else if (type.first == info_t::zero_one)
		decode_info_one(v_compressed, v_data);
	else if (type.first == info_t::constant)
		decode_info_constant(type.second, v_compressed, v_data);
	else if (type.first == info_t::zero_constant)
		decode_info_constant(type.second, v_compressed, v_data);
	else
		decode_info_any(v_size, v_compressed, v_data);
}

// EOF
