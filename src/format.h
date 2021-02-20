#pragma once
// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include <cstdint>
#include <array>
#include <vector>
#include <unordered_map>
#include "utils.h"
#include "io.h"
#include "rc.h"
#include "sub_rc.h"
#include "context_hm.h"
#include "hm.h"

using namespace std;

// ************************************************************************************
class CFormatCompress
{
	enum class info_t {unknown, zero, one, zero_one, constant, zero_constant, any};

	string desc;
	uint32_t no_samples;
	uint32_t compression_level;

	CVectorIOStream* vios_i;
	CVectorIOStream* vios_o;

	vector<uint8_t> v_vios_i;
	vector<uint8_t> v_vios_o;

	CRangeEncoder<CVectorIOStream>* rce;
	CRangeDecoder<CVectorIOStream>* rcd;

	using ModelType_16_19_16 = CSimpleModel<16, 19, 16>;
	using ModelType_2_15_1 = CSimpleModel<2, 15, 1>;
	using ModelType_256_16_1 = CAdjustableModelEmb<256, 16, 1>;
	using ModelType_256_19_128 = CAdjustableModelEmb<256, 19, 128>;
	using ModelType_2_19_16 = CSimpleModel<2, 19, 16>;

	using ctx_map_16_19_16_t = CContextHM<CRangeCoderModel<ModelType_16_19_16, CVectorIOStream, 16, 19, 16>>;
	using ctx_map_2_15_1_t = CContextHM<CRangeCoderModel<ModelType_2_15_1, CVectorIOStream, 2, 15, 1>>;
	using ctx_map_256_16_1_t = CContextHM<CRangeCoderModel<ModelType_256_16_1, CVectorIOStream, 256, 16, 1>>;
	using ctx_map_256_19_128_t = CContextHM<CRangeCoderModel<ModelType_256_19_128, CVectorIOStream, 256, 19, 128>>;
	using ctx_map_2_19_16_t = CContextHM<CRangeCoderModel<ModelType_2_19_16, CVectorIOStream, 2, 19, 16>>;

	ctx_map_2_19_16_t ctx_map_same;
	ctx_map_2_15_1_t ctx_map_known;
	ctx_map_2_19_16_t ctx_map_known2;
	ctx_map_256_16_1_t ctx_map_plain;
	ctx_map_256_19_128_t ctx_map_code;
	ctx_map_16_19_16_t ctx_map_entropy_type;

	pair<info_t, uint32_t> type = { info_t::unknown, 0 };
	uint32_t ctx_mode = 0;
	uint32_t max_dict_size = 1u << 20;
	const uint32_t ht_empty_key = 0x7fffffffu;
	hash_map_lp<uint32_t, uint32_t, std::equal_to<uint32_t>, MurMur32Hash> dict;
	vector<uint32_t> dict_dec;

	pair<info_t, uint32_t> determine_info_type(vector<uint32_t>& v_size);

	template<unsigned NO_SYMBOLS, unsigned MAX_LOG_COUNTER, unsigned ADDER>
	CRangeCoderModel<CSimpleModel<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>* 
		find_rce_coder(CContextHM<CRangeCoderModel<CSimpleModel<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>>& map, context_t ctx)
	{
		auto p = map.find(ctx);

		if (p == nullptr)
			map.insert(ctx, p = new CRangeCoderModel<CSimpleModel<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>(rce, nullptr, true));

		return p;
	}

	template<unsigned NO_SYMBOLS, unsigned MAX_LOG_COUNTER, unsigned ADDER>
	CRangeCoderModel<CAdjustableModel<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>* 
		find_rce_coder(CContextHM<CRangeCoderModel<CAdjustableModel<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>>& map, context_t ctx)
	{
		auto p = map.find(ctx);

		if (p == nullptr)
			map.insert(ctx, p = new CRangeCoderModel<CAdjustableModel<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>(rce, nullptr, true));

		return p;
	}

	template<unsigned NO_SYMBOLS, unsigned MAX_LOG_COUNTER, unsigned ADDER>
	CRangeCoderModel<CAdjustableModelEmb<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>* 
		find_rce_coder(CContextHM<CRangeCoderModel<CAdjustableModelEmb<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>>& map, context_t ctx)
	{
		auto p = map.find(ctx);

		if (p == nullptr)
			map.insert(ctx, p = new CRangeCoderModel<CAdjustableModelEmb<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>(rce, nullptr, true));

		return p;
	}

	template<unsigned NO_SYMBOLS, unsigned MAX_LOG_COUNTER, unsigned ADDER>
	CRangeCoderModel<CSimpleModel<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>*
		find_rcd_coder(CContextHM<CRangeCoderModel<CSimpleModel<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>>& map, context_t ctx)
	{
		auto p = map.find(ctx);

		if (p == nullptr)
			map.insert(ctx, p = new CRangeCoderModel<CSimpleModel<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>(rcd, nullptr, false));

		return p;
	}

	template<unsigned NO_SYMBOLS, unsigned MAX_LOG_COUNTER, unsigned ADDER>
	CRangeCoderModel<CAdjustableModel<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>*
		find_rcd_coder(CContextHM<CRangeCoderModel<CAdjustableModel<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>>& map, context_t ctx)
	{
		auto p = map.find(ctx);

		if (p == nullptr)
			map.insert(ctx, p = new CRangeCoderModel<CAdjustableModel<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>(rcd, nullptr, false));

		return p;
	}

	template<unsigned NO_SYMBOLS, unsigned MAX_LOG_COUNTER, unsigned ADDER>
	CRangeCoderModel<CAdjustableModelEmb<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>*
		find_rcd_coder(CContextHM<CRangeCoderModel<CAdjustableModelEmb<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>>& map, context_t ctx)
	{
		auto p = map.find(ctx);

		if (p == nullptr)
			map.insert(ctx, p = new CRangeCoderModel<CAdjustableModelEmb<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>(rcd, nullptr, false));

		return p;
	}


	template<unsigned NO_SYMBOLS, unsigned MAX_LOG_COUNTER, unsigned ADDER>
	void prefetch_rc_coder(CContextHM<CRangeCoderModel<CAdjustableModel<NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>, CVectorIOStream, NO_SYMBOLS, MAX_LOG_COUNTER, ADDER>>& map, context_t ctx)
	{
		map.prefetch(ctx);
	}



	template <unsigned SIZE> 
	class array_hash
	{
	public:
		size_t operator()(const array<uint32_t, SIZE>& x) const {
			size_t r = 0;
			for (auto a : x)
				r ^= hash<uint32_t>{}(a);

			return r;
		}
	};

	void encode_ctx_type(int ctx_type);
	int decode_ctx_type();

	template <unsigned SIZE> double entropy(vector<array<uint32_t, SIZE>> &vec);
	double entropy_0(vector<pair<uint32_t, uint32_t>>& vec)
	{
		double r = (double)vec.size() * max(8.0, log2(vec.size()));

		double sum = 0.0;

		for (auto& x : vec)
			sum += x.second;
		for (auto& x : vec)
			r -= (double)x.second * log2((double)x.second / sum);

		return r;
	}

	void encode_info_one(vector<uint8_t>& v_data, vector<uint8_t>& v_compressed);
	void encode_info_constant(uint32_t val, vector<uint8_t>& v_data, vector<uint8_t>& v_compressed);
	void encode_info_any(vector<uint32_t>& v_size, vector<uint8_t>& v_data, vector<uint8_t>& v_compressed);

	void decode_info_one(vector<uint8_t>& v_compressed, vector<uint8_t>& v_data);
	void decode_info_constant(uint32_t val, vector<uint8_t>& v_compressed, vector<uint8_t>& v_data);
	void decode_info_any(vector<uint32_t>& v_size, vector<uint8_t>& v_compressed, vector<uint8_t>& v_data);

	void encode_plain(uint32_t x);
	void encode_code(context_t ctx, uint32_t code);

	void update_code_enc(context_t ctx, uint32_t code);
	void update_code_dec(context_t ctx, uint32_t code);

	uint32_t decode_plain();
	uint32_t decode_code(context_t ctx);

	void encode_format_one(vector<uint32_t>& v_size, vector<uint8_t>& v_data, vector<uint8_t>& v_compressed);
	void encode_format_many(vector<uint32_t>& v_size, vector<uint8_t>& v_data, vector<uint8_t>& v_compressed);

	void decode_format_one(vector<uint32_t>& v_size, vector<uint8_t>& v_data, vector<uint8_t>& v_compressed);
	void decode_format_many(vector<uint32_t>& v_size, vector<uint8_t>& v_data, vector<uint8_t>& v_compressed);

	void append_uint32(vector<uint8_t>& vec, uint32_t x)
	{
		vec.emplace_back(x & 0xff);
		vec.emplace_back((x >> 8) & 0xff);
		vec.emplace_back((x >> 16) & 0xff);
		vec.emplace_back(x >> 24);
	}

public:
	CFormatCompress(string _desc, uint32_t _compression_level);
	~CFormatCompress();

	void SetNoSamples(uint32_t _no_samples);

	void EncodeFormat(vector<uint32_t>& v_size, vector<uint8_t>& v_data, vector<uint8_t>& v_compressed);
	void EncodeInfo(vector<uint32_t>& v_size, vector<uint8_t>& v_data, vector<uint8_t>& v_compressed);
	
	void DecodeFormat(vector<uint32_t>& v_size, vector<uint8_t>& v_compressed, vector<uint8_t>& v_data);
	void DecodeInfo(vector<uint32_t>& v_size, vector<uint8_t>& v_compressed, vector<uint8_t>& v_data);
};

// EOF
