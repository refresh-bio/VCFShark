#pragma once
// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2020-12-18
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

using namespace std;

// ************************************************************************************
class CFormatCompress
{
	enum class info_t {unknown, zero, one, zero_one, constant, zero_constant, any};

	uint32_t no_samples;

	CVectorIOStream* vios_i;
	CVectorIOStream* vios_o;

	vector<uint8_t> v_vios_i;
	vector<uint8_t> v_vios_o;

	CRangeEncoder<CVectorIOStream>* rce;
	CRangeDecoder<CVectorIOStream>* rcd;

	using ModelType = CAdjustableModel;

	typedef CContextHM<CRangeCoderModel<ModelType, CVectorIOStream>> ctx_map_t;

	ctx_map_t ctx_map_same;
	ctx_map_t ctx_map_known;
	ctx_map_t ctx_map_plain;
	ctx_map_t ctx_map_code;
	ctx_map_t ctx_map_entropy_type;

	pair<info_t, uint32_t> type = { info_t::unknown, 0 };
	uint32_t ctx_mode = 0;
	uint32_t max_dict_size = 1u << 20;
	unordered_map<uint32_t, uint32_t> dict;

	pair<info_t, uint32_t> determine_info_type(vector<uint32_t>& v_size);

	inline ctx_map_t::value_type find_rce_coder(ctx_map_t &map, context_t ctx, uint32_t no_symbols, uint32_t max_log_counter, uint32_t adder);
	inline ctx_map_t::value_type find_rcd_coder(ctx_map_t& map, context_t ctx, uint32_t no_symbols, uint32_t max_log_counter, uint32_t adder);

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

	void append_uint32(vector<uint8_t>& vec, uint32_t x);

public:
	CFormatCompress();
	~CFormatCompress();

	void SetNoSamples(uint32_t _no_samples);

	void EncodeFormat(vector<uint32_t>& v_size, vector<uint8_t>& v_data, vector<uint8_t>& v_compressed);
	void EncodeInfo(vector<uint32_t>& v_size, vector<uint8_t>& v_data, vector<uint8_t>& v_compressed);
	
	void DecodeFormat(vector<uint32_t>& v_size, vector<uint8_t>& v_compressed, vector<uint8_t>& v_data);
	void DecodeInfo(vector<uint32_t>& v_size, vector<uint8_t>& v_compressed, vector<uint8_t>& v_data);
};

// EOF
