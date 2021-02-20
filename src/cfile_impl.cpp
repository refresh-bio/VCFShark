// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include <memory>
#include <iostream>
#include <set>
#include <future>
using namespace std;

#include "cfile.h"
#include "utils.h"

// ************************************************************************************
bool CCompressedFile::load_descriptions()
{
	string s;

	// Load file header and characteristics
	vector<uint8_t> v_desc;
	size_t p_desc = 0;
	uint64_t tmp;

	auto stream_id = archive->GetStreamId("db_params");
	if (stream_id < 0)
	{
		std::cerr << "Corrupted archive!\n";
		exit(0);
	}

	size_t aux;

	archive->GetPart(stream_id, v_desc, aux);
	read(v_desc, p_desc, no_variants);
	read(v_desc, p_desc, no_samples);
	read_fixed(v_desc, p_desc, tmp, 1);
	ploidy = (uint8_t)tmp;
	read(v_desc, p_desc, neglect_limit);
	i_variant = 0;

	read(v_desc, p_desc, no_keys);

	uint32_t tmp32;
	read(v_desc, p_desc, tmp32);
	gt_key_id = (int)tmp32;

	keys.resize(no_keys);

	for (uint32_t i = 0; i < no_keys; ++i)
	{
		read(v_desc, p_desc, keys[i].key_id);
		read_fixed(v_desc, p_desc, tmp, 1);
		keys[i].keys_type = (key_type_t)tmp;
		read_fixed(v_desc, p_desc, tmp, 1);
		keys[i].type = (int8_t) tmp;
	}

	// Load variant descriptions
	for (auto d : {
		make_tuple(ref(v_rd_meta), ref(v_cd_meta), ref(p_meta), 4, "meta"),
		make_tuple(ref(v_rd_header), ref(v_cd_header), ref(p_header), 4, "header"),
		make_tuple(ref(v_rd_samples), ref(v_cd_samples), ref(p_samples), 4, "samples")
		})
	{
		stream_id = archive->GetStreamId("db_" + string(get<4>(d)));

		archive->GetPart(stream_id, get<1>(d), aux);

		CBSCWrapper bsc;
		vector<uint8_t> v_tmp;

		bsc.InitDecompress();
		bsc.Decompress(get<1>(d), v_tmp);

		get<0>(d) = move(v_tmp);
		get<2>(d) = 0;
	}

	v_meta.clear();
	read(v_rd_meta, p_meta, v_meta);

	v_header.clear();
	read(v_rd_header, p_header, v_header);

	v_samples.clear();
	string sample;
	for (uint32_t i = 0; i < no_samples; ++i)
	{
		read(v_rd_samples, p_samples, sample);
		v_samples.emplace_back(sample);
	}

	return true;
}

// ************************************************************************************
bool CCompressedFile::save_descriptions()
{
	vector<uint8_t> v_desc;

	append(v_desc, no_variants);
	append(v_desc, no_samples);
	append_fixed(v_desc, ploidy, 1);
	append(v_desc, neglect_limit);

	// Save key descriptions
	append(v_desc, no_keys);
	append(v_desc, (int)gt_key_id);

	for (uint32_t i = 0; i < no_keys; ++i)
	{
		append(v_desc, keys[i].key_id);
		append_fixed(v_desc, static_cast<uint64_t>(keys[i].keys_type), 1);
		append_fixed(v_desc, keys[i].type, 1);
	}

	auto stream_id = archive->RegisterStream("db_params");
	archive->AddPart(stream_id, v_desc);
	archive->SetRawSize(stream_id, v_desc.size());

	append(v_rd_meta, v_meta);
	append(v_rd_header, v_header);

	for (auto& x : v_samples)
		append(v_rd_samples, x);

	// Save variant descriptions
	for (auto d : {
		make_tuple(ref(v_rd_meta), ref(v_cd_meta), 4, "meta"),
		make_tuple(ref(v_rd_header), ref(v_cd_header), 4, "header"),
		make_tuple(ref(v_rd_samples), ref(v_cd_samples), 4, "samples")
		})
	{
		CBSCWrapper bsc;
		size_t r_size = get<0>(d).size();
		get<1>(d).clear();

		bsc.InitCompress(p_bsc_meta);
		bsc.Compress(get<0>(d), get<1>(d));

		auto stream_id = archive->RegisterStream("db_" + string(get<3>(d)));
		archive->AddPart(stream_id, get<1>(d));
		archive->SetRawSize(stream_id, r_size);
	}

	return true;
}

// ************************************************************************************
void CCompressedFile::lock_coder_compressor(SPackage& pck)
{
	unique_lock<mutex> lck(mtx_v_coder);
	cv_v_coder.wait(lck, [&, this] {
		int sid = pck.key_id;
		if (pck.type == SPackage::package_t::db)
			sid = no_keys + pck.db_id;

		return (int) v_coder_part_ids[sid] == pck.part_id;
		});
}

// ************************************************************************************
bool CCompressedFile::check_coder_compressor(SPackage& pck)
{
	unique_lock<mutex> lck(mtx_v_coder);
	int sid = pck.key_id;
	if (pck.type == SPackage::package_t::db)
		sid = no_keys + pck.db_id;

	return (int) v_coder_part_ids[sid] == pck.part_id;
}

// ************************************************************************************
void CCompressedFile::unlock_coder_compressor(SPackage& pck)
{
	lock_guard<mutex> lck(mtx_v_coder);
	int sid = pck.key_id;
	if (pck.type == SPackage::package_t::db)
		sid = no_keys + pck.db_id;

	++v_coder_part_ids[sid];
	cv_v_coder.notify_all();
}

// ************************************************************************************
void CCompressedFile::lock_text_compressor(SPackage& pck)
{
	unique_lock<mutex> lck(mtx_v_text);
	cv_v_text.wait(lck, [&, this] {
		int sid = pck.key_id;
		if (pck.type == SPackage::package_t::db)
			sid = no_keys + pck.db_id;

		return (int) v_text_part_ids[sid] == pck.part_id;
		});
}

// ************************************************************************************
void CCompressedFile::unlock_text_compressor(SPackage& pck)
{
	lock_guard<mutex> lck(mtx_v_text);
	int sid = pck.key_id;
	if (pck.type == SPackage::package_t::db)
		sid = no_keys + pck.db_id;

	++v_text_part_ids[sid];
	cv_v_text.notify_all();
}

// ************************************************************************************
void CCompressedFile::skip_text_compressor(SPackage& pck)
{
	unique_lock<mutex> lck(mtx_v_text);
	int sid = pck.key_id;
	if (pck.type == SPackage::package_t::db)
		sid = no_keys + pck.db_id;

	cv_v_text.wait(lck, [&, this] {
		return (int) v_text_part_ids[sid] == pck.part_id;
		});

	++v_text_part_ids[sid];
	cv_v_text.notify_all();
}

// ************************************************************************************
void CCompressedFile::compress_field(SPackage& pck, vector<uint8_t>& v_compressed, vector<uint8_t>& v_tmp)
{
	CBSCWrapper* bsc_size = v_bsc_size[pck.key_id];
	CBSCWrapper* bsc_data = v_bsc_data[pck.key_id];
	size_t raw_size;

	if (pck.v_data.size())
	{
		bool is_pp_compressed = false;

		if (keys[pck.key_id].type == BCF_HT_STR && 64 * pck.v_size.size() < pck.v_data.size())
		{
			vector<uint8_t> v_pp;
			lock_text_compressor(pck);
			v_text_pp[pck.key_id].EncodeText(pck.v_data, v_pp);
			unlock_text_compressor(pck);
			lock_coder_compressor(pck);
			bsc_data->Compress(v_pp, v_compressed);
			raw_size = v_pp.size();

			is_pp_compressed = true;
		}
		else
		{
			skip_text_compressor(pck);
			lock_coder_compressor(pck);
			bsc_data->Compress(pck.v_data, v_compressed);
			raw_size = pck.v_data.size();
		}

		archive->AddPartComplete(pck.stream_id_data, pck.part_id, v_compressed, raw_size + (is_pp_compressed ? pp_compress_flag : 0));
	}
	else
	{
		v_compressed.clear();
		skip_text_compressor(pck);
		archive->AddPartComplete(pck.stream_id_data, pck.part_id, v_compressed, pck.v_data.size());
	}

	v_tmp.resize(pck.v_size.size() * 4);
	copy_n((uint8_t*)pck.v_size.data(), v_tmp.size(), v_tmp.data());

	bsc_size->Compress(v_tmp, v_compressed);
	archive->AddPartComplete(pck.stream_id_size, pck.part_id, v_compressed, pck.v_size.size());

	unlock_coder_compressor(pck);
}

// ************************************************************************************
void CCompressedFile::decompress_field(SPackage* pck, size_t raw_size, vector<uint8_t>& v_tmp)
{
	if (pck->is_func)
	{
		load_function("func_" + to_string(pck->key_id) + "_data", pck->stream_id_src, pck->fun);
		v_packages[pck->key_id] = pck;

		return;
	}

	pck->stream_id_data = archive->GetStreamId("key_" + to_string(pck->key_id) + "_data");

	CBSCWrapper* bsc_size = v_bsc_size[pck->key_id];
	CBSCWrapper* bsc_data = v_bsc_data[pck->key_id];

	bsc_size->Decompress(pck->v_compressed, v_tmp);

	pck->v_size.resize(raw_size);
	copy_n(v_tmp.data(), raw_size * 4, (uint8_t*)pck->v_size.data());

	archive->GetPart(pck->stream_id_data, pck->v_compressed, raw_size);

	bool is_pp_compressed = false;

	if (raw_size >= pp_compress_flag)
	{
		raw_size -= pp_compress_flag;
		is_pp_compressed = true;
	}

	pck->v_data.resize(raw_size);

	if (raw_size)
	{
		bsc_data->Decompress(pck->v_compressed, pck->v_data);

		if (is_pp_compressed)
		{
			vector<uint8_t> v_decompressed;
			v_text_pp[pck->key_id].DecodeText(pck->v_data, v_decompressed);
			swap(pck->v_data, v_decompressed);
		}
	}
}

// ************************************************************************************
void CCompressedFile::compress_format(SPackage& pck, vector<uint8_t>& v_compressed, vector<uint8_t>& v_tmp)
{
	CBSCWrapper* bsc_size = v_bsc_size[pck.key_id];
	CFormatCompress* format_compress = v_format_compress[pck.key_id];
//	size_t raw_size;

	if (pck.v_data.size())
	{
		skip_text_compressor(pck);
		lock_coder_compressor(pck);
		format_compress->EncodeFormat(pck.v_size, pck.v_data, v_compressed);

		archive->AddPartComplete(pck.stream_id_data, pck.part_id, v_compressed, pck.v_data.size());
	}
	else
	{
		v_compressed.clear();
		skip_text_compressor(pck);
		archive->AddPartComplete(pck.stream_id_data, pck.part_id, v_compressed, pck.v_data.size());
	}

	v_tmp.resize(pck.v_size.size() * 4);
	copy_n((uint8_t*)pck.v_size.data(), v_tmp.size(), v_tmp.data());

	bsc_size->Compress(v_tmp, v_compressed);
	archive->AddPartComplete(pck.stream_id_size, pck.part_id, v_compressed, pck.v_size.size());

	unlock_coder_compressor(pck);
}

// ************************************************************************************
void CCompressedFile::decompress_format(SPackage* pck, size_t raw_size, vector<uint8_t>& v_tmp)
{
	if (pck->is_func)
	{
		load_function("func_" + to_string(pck->key_id) + "_data", pck->stream_id_src, pck->fun);
		v_packages[pck->key_id] = pck;

		return;
	}

	pck->stream_id_data = archive->GetStreamId("key_" + to_string(pck->key_id) + "_data");

	CBSCWrapper* bsc_size = v_bsc_size[pck->key_id];
	CFormatCompress* format_compress = v_format_compress[pck->key_id];

	bsc_size->Decompress(pck->v_compressed, v_tmp);

	pck->v_size.resize(raw_size);
	copy_n(v_tmp.data(), raw_size * 4, (uint8_t*)pck->v_size.data());

	archive->GetPart(pck->stream_id_data, pck->v_compressed, raw_size);

//	bool is_pp_compressed = false;

	if (raw_size >= pp_compress_flag)
	{
		raw_size -= pp_compress_flag;
//		is_pp_compressed = true;
	}

	pck->v_data.resize(raw_size);

	if (raw_size)
		format_compress->DecodeFormat(pck->v_size, pck->v_compressed, pck->v_data);
}

// ************************************************************************************
void CCompressedFile::compress_info(SPackage& pck, vector<uint8_t>& v_compressed, vector<uint8_t>& v_tmp)
{
	CBSCWrapper* bsc_size = v_bsc_size[pck.key_id];
	CFormatCompress* format_compress = v_format_compress[pck.key_id];
//	size_t raw_size;

	if (pck.v_data.size())
	{
		skip_text_compressor(pck);
		lock_coder_compressor(pck);
		format_compress->EncodeInfo(pck.v_size, pck.v_data, v_compressed);

		archive->AddPartComplete(pck.stream_id_data, pck.part_id, v_compressed, pck.v_data.size());
	}
	else
	{
		v_compressed.clear();
		skip_text_compressor(pck);
		archive->AddPartComplete(pck.stream_id_data, pck.part_id, v_compressed, pck.v_data.size());
	}

	v_tmp.resize(pck.v_size.size() * 4);
	copy_n((uint8_t*)pck.v_size.data(), v_tmp.size(), v_tmp.data());

	bsc_size->Compress(v_tmp, v_compressed);
	archive->AddPartComplete(pck.stream_id_size, pck.part_id, v_compressed, pck.v_size.size());

	unlock_coder_compressor(pck);
}

// ************************************************************************************
void CCompressedFile::decompress_info(SPackage* pck, size_t raw_size, vector<uint8_t>& v_tmp)
{
	if (pck->is_func)
	{
		load_function("func_" + to_string(pck->key_id) + "_data", pck->stream_id_src, pck->fun);
		v_packages[pck->key_id] = pck;

		return;
	}

	pck->stream_id_data = archive->GetStreamId("key_" + to_string(pck->key_id) + "_data");

	CBSCWrapper* bsc_size = v_bsc_size[pck->key_id];
	CFormatCompress* format_compress = v_format_compress[pck->key_id];

	bsc_size->Decompress(pck->v_compressed, v_tmp);

	pck->v_size.resize(raw_size);
	copy_n(v_tmp.data(), raw_size * 4, (uint8_t*)pck->v_size.data());

	archive->GetPart(pck->stream_id_data, pck->v_compressed, raw_size);

	pck->v_data.resize(raw_size);

	if (raw_size)
		format_compress->DecodeInfo(pck->v_size, pck->v_compressed, pck->v_data);
}

// ************************************************************************************
void CCompressedFile::compress_db(SPackage& pck, vector<uint8_t>& v_compressed, vector<uint8_t>& v_tmp)
{
	CBSCWrapper* bsc_size = v_bsc_db_size[pck.db_id];
	CBSCWrapper* bsc_data = v_bsc_db_data[pck.db_id];
	size_t raw_size;

	v_tmp.resize(pck.v_size.size() * 4);
	copy_n((uint8_t*)pck.v_size.data(), v_tmp.size(), v_tmp.data());

	lock_coder_compressor(pck);

	bsc_size->Compress(v_tmp, v_compressed);
	archive->AddPartComplete(pck.stream_id_size, pck.part_id, v_compressed, pck.v_size.size());

	if (pck.v_data.size())
	{
		raw_size = pck.v_data.size();
		bsc_data->Compress(pck.v_data, v_compressed);

		archive->AddPartComplete(pck.stream_id_data, pck.part_id, v_compressed, raw_size);
	}
	else
	{
		v_compressed.clear();
		archive->AddPartComplete(pck.stream_id_data, pck.part_id, v_compressed, pck.v_data.size());
	}

	unlock_coder_compressor(pck);
}

// ************************************************************************************
void CCompressedFile::decompress_db(SPackage* pck, size_t raw_size, vector<uint8_t>& v_tmp)
{
	CBSCWrapper* bsc_size = v_bsc_db_size[pck->db_id];
	CBSCWrapper* bsc_data = v_bsc_db_data[pck->db_id];

	bsc_size->Decompress(pck->v_compressed, v_tmp);

	pck->v_size.resize(raw_size);
	copy_n(v_tmp.data(), raw_size * 4, (uint8_t*)pck->v_size.data());

	archive->GetPart(pck->stream_id_data, pck->v_compressed, raw_size);

	pck->v_data.resize(raw_size);

	if (raw_size)
		bsc_data->Decompress(pck->v_compressed, pck->v_data);
}

#if 1
// ************************************************************************************
void CCompressedFile::compress_gt(SPackage& pck)
{
	int i_vec = 0;
	vector<uint32_t> v_tmp_reo;
	vector<pair<uint32_t, uint32_t>> v_rle;

	vector<uint32_t> v_res;

	lock_coder_compressor(pck);

	// *** Reorganization of haplotypes
	for (size_t i = 0; i < pck.v_data.size(); i += pck.v_size[i_vec++] * 4)
	{
		uint32_t* vec = (uint32_t*)(pck.v_data.data() + i);
		uint32_t no_haplotypes = pck.v_size[i_vec] / no_samples;
		uint32_t max_gt_val = 0;

		v_tmp_reo.resize(pck.v_size[i_vec]);

		// Change of status of the 1st haplotype
		if (no_haplotypes > 1)
			for (uint32_t k = 0; k < no_samples; ++k)
				if (vec[k * no_haplotypes + 1] & 1)
					vec[k * no_haplotypes] += 1;

		for (uint32_t j = 0; j < no_haplotypes; ++j)
			for (uint32_t k = 0; k < no_samples; ++k)
			{
				uint32_t gt_val = vec[k * no_haplotypes + j];
				if (gt_val == 0x80000001u)
					gt_val = 0;
				else
					++gt_val;

				v_tmp_reo[j * no_samples + k] = gt_val;

				if (gt_val > max_gt_val)
				{
					max_gt_val = gt_val;
#ifdef LOG_INFO
					if (max_gt_val > 64)
						cout << "***\n";
#endif
				}
			}

		pbwt.EncodeFlexible(max_gt_val, v_tmp_reo, v_rle);

		v_rle.back().second = 0;

		for (auto x : v_rle)
		{
			v_res.emplace_back(x.first);
			v_res.emplace_back(x.second);
		}
	}

	if (pck.v_data.size())
	{
		for (auto& x : pck.v_size)
			x /= no_samples;

		//  BSC compression
		CBSCWrapper* bsc_size = v_bsc_size[pck.key_id];

		vector<uint8_t> v_tmp;
		vector<uint8_t> v_compressed;
		size_t raw_size;

		v_tmp.resize(pck.v_size.size() * 4);
		copy_n((uint8_t*)pck.v_size.data(), v_tmp.size(), v_tmp.data());

		bsc_size->Compress(v_tmp, v_compressed);
		archive->AddPartComplete(pck.stream_id_size, pck.part_id, v_compressed, pck.v_size.size());

		// RC compression
		v_vios_o.clear();
		rce->Start();
		ctx_prefix = context_prefix_mask;
		ctx_symbol = context_symbol_mask;

		for (size_t i = 0; i < v_res.size(); i += 2)
		{
			encode_run_len(v_res[i], v_res[i + 1]);
			if (v_res[i + 1] == 0)
			{
				ctx_prefix = context_prefix_mask;
				ctx_symbol = context_symbol_mask;
			}
		}

		rce->End();

		raw_size = v_res.size();

		archive->AddPartComplete(pck.stream_id_data, pck.part_id, v_vios_o, raw_size);
	}
	else
	{
		vector<uint8_t> v_compressed;
		archive->AddPartComplete(pck.stream_id_data, pck.part_id, v_compressed, 0);
	}

	unlock_coder_compressor(pck);
}
#endif

#if 0
// ************************************************************************************
void CCompressedFile::compress_gt(SPackage& pck)
{
	int i_vec = 0;
	vector<uint32_t> v_tmp_reo;
	vector<pair<uint32_t, uint32_t>> v_rle;

	vector<uint32_t> v_res;

	if (pck.v_data.empty())
	{
		vector<uint8_t> v_compressed;
		archive->AddPartComplete(pck.stream_id_data, pck.part_id, v_compressed, 0);

		return;
	}

	lock_coder_compressor(pck);

	// *** Reorganization of haplotypes
	for (int i = 0; i < pck.v_data.size(); i += pck.v_size[i_vec++] * 4)
	{
		uint32_t* vec = (uint32_t*)(pck.v_data.data() + i);
		int no_haplotypes = pck.v_size[i_vec] / no_samples;
		uint32_t max_gt_val = 0;

		v_tmp_reo.resize(pck.v_size[i_vec]);

		// Change of status of the 1st haplotype
		if (no_haplotypes > 1)
			for (int k = 0; k < no_samples; ++k)
				if (vec[k * no_haplotypes + 1] & 1)
					vec[k * no_haplotypes] += 1;

		for (int j = 0; j < no_haplotypes; ++j)
			for (int k = 0; k < no_samples; ++k)
			{
				uint32_t gt_val = vec[k * no_haplotypes + j];
				if (gt_val == 0x80000001u)
					gt_val = 0;
				else
					++gt_val;

				v_tmp_reo[j * no_samples + k] = gt_val;

				if (gt_val > max_gt_val)
				{
					max_gt_val = gt_val;
#ifdef LOG_INFO
					if (max_gt_val > 64)
						cout << "***\n";
#endif
				}
			}

		pbwt.EncodeFlexible(max_gt_val, v_tmp_reo, v_rle);

		v_rle.back().second = 0;

		for (auto x : v_rle)
		{
			v_res.emplace_back(x.first);
			v_res.emplace_back(x.second);
		}
	}

	for (auto& x : pck.v_size)
		x /= no_samples;

	size_t raw_size;

	//  BSC compression
	auto bsc_async = async([&] {
		CBSCWrapper* bsc_size = v_bsc_size[pck.key_id];

		vector<uint8_t> v_tmp;
		vector<uint8_t> v_compressed;

		v_tmp.resize(pck.v_size.size() * 4);
		copy_n((uint8_t*)pck.v_size.data(), v_tmp.size(), v_tmp.data());

		bsc_size->Compress(v_tmp, v_compressed);
		archive->AddPartComplete(pck.stream_id_size, pck.part_id, v_compressed, pck.v_size.size());
		});

	// RC compression
	v_vios_o.clear();
	rce->Start();
	ctx_prefix = context_prefix_mask;
	ctx_symbol = context_symbol_mask;

	for (int i = 0; i < v_res.size(); i += 2)
	{
		encode_run_len(v_res[i], v_res[i + 1]);
		if (v_res[i + 1] == 0)
		{
			ctx_prefix = context_prefix_mask;
			ctx_symbol = context_symbol_mask;
		}
	}

	rce->End();

	bsc_async.wait();

	raw_size = v_res.size();

	archive->AddPartComplete(pck.stream_id_data, pck.part_id, v_vios_o, raw_size);

	unlock_coder_compressor(pck);
}
#endif

// ************************************************************************************
void CCompressedFile::decompress_gt(SPackage* pck, size_t raw_size)
{
	CBSCWrapper* bsc_size = v_bsc_size[pck->key_id];
//	CBSCWrapper* bsc_data = v_bsc_data[pck->key_id];

	if (raw_size == 0)
	{
		pck->v_data.clear();
		pck->v_size.clear();

		return;
	}

	vector<uint8_t> v_tmp;
	bsc_size->Decompress(pck->v_compressed, v_tmp);

	pck->v_size.resize(raw_size);
	copy_n(v_tmp.data(), raw_size * 4, (uint8_t*)pck->v_size.data());

	pck->stream_id_data = archive->GetStreamId("key_" + to_string(pck->key_id) + "_data");

	archive->GetPart(pck->stream_id_data, pck->v_compressed, raw_size);
	pck->v_data.resize(raw_size);

	vector<pair<uint32_t, uint32_t>> v_full_rle;

	if (raw_size)
	{
		v_vios_i = pck->v_compressed;
		vios_i->RestartRead();

		rcd->Start();
		ctx_prefix = context_prefix_mask;
		ctx_symbol = context_symbol_mask;

		size_t cur_size = 0;
		int i_variant = 0;

		uint32_t symbol;
		uint32_t len;
		uint32_t cur_variant_size = 0;

		for (size_t i = 0; i < raw_size; i += 2)
		{
			decode_run_len(symbol, len);
			if (len == 0)
			{
				ctx_prefix = context_prefix_mask;
				ctx_symbol = context_symbol_mask;
				len = pck->v_size[i_variant++] * no_samples - cur_variant_size;
				cur_variant_size = 0;
			}
			else
				cur_variant_size += len;

			cur_size += len;

			v_full_rle.emplace_back(symbol, len);
		}

		rcd->End();
	}

	// PBWT decoding
	size_t total_data_size = 0;
	for (auto& x : pck->v_size)
	{
		x *= no_samples;

		total_data_size += x;
	}

	pck->v_data.resize(total_data_size * 4);

	size_t i_data = 0;

	int i_variant = 0;
	vector<pair<uint32_t, uint32_t>> v_rle;
	vector<uint32_t> v_output, vec;

	for (size_t i = 0; i < v_full_rle.size();)
	{
		uint32_t c_variant_len = 0;
		uint32_t variant_size = pck->v_size[i_variant++];

		v_rle.clear();

		uint32_t max_val = 0;

		for (; c_variant_len < variant_size; ++i)
		{
			v_rle.emplace_back(v_full_rle[i]);
			c_variant_len += v_full_rle[i].second;

			if (v_full_rle[i].first > max_val)
				max_val = v_full_rle[i].first;
		}

		pbwt.DecodeFlexible(max_val, v_rle, v_output);

		vec.resize(v_output.size());

		uint32_t no_haplotypes = (uint32_t) (v_output.size() / no_samples);

		for (uint32_t j = 0; j < no_haplotypes; ++j)
			for (uint32_t k = 0; k < no_samples; ++k)
			{
				uint32_t gt_val = v_output[j * no_samples + k];

				if (gt_val == 0)
					gt_val = 0x80000001u;
				else
					--gt_val;

				vec[k * no_haplotypes + j] = gt_val;
			}

		// Recovery of the 1st haplotype status
		if (no_haplotypes > 1)
			for (uint32_t k = 0; k < no_samples; ++k)
				if (vec[k * no_haplotypes] & 1)
					vec[k * no_haplotypes] -= 1;


		copy_n((uint8_t*)vec.data(), vec.size() * 4, pck->v_data.data() + i_data);
		i_data += v_output.size() * 4;
	}
}

// ************************************************************************************
void CCompressedFile::encode_run_len(uint32_t symbol, uint32_t len)
{
	// Encode symbol
	auto rc_sym = find_rce_coder(rcd_coders_rl_sym, ctx_symbol + context_symbol_flag);

	if (symbol < 15)
		rc_sym->Encode(symbol);
	else
	{
		uint32_t x = symbol;

		while (true)
		{
			if (x < 15)
			{
				rc_sym->Encode(x);
				break;
			}

			rc_sym->Encode(15);
			x -= 15;
		}

		symbol = 15;
	}

	ctx_symbol <<= 4;
	ctx_symbol += symbol;
	ctx_symbol &= context_symbol_mask;

	rce_coders_rl_sym.prefetch(ctx_symbol + context_symbol_flag);

	ctx_prefix <<= 4;
	ctx_prefix += (context_t)symbol;
	ctx_prefix &= context_prefix_mask;

	// Encode run length
	auto rc_p = find_rce_coder(rce_coders_rl_pref, ctx_prefix + context_prefix_flag);

	uint32_t prefix = ilog2(len);

	ctx_prefix <<= 4;
	ctx_prefix += (context_t)prefix;
	ctx_prefix &= context_prefix_mask;

	rce_coders_rl_pref.prefetch(ctx_prefix + context_prefix_flag);

	if (prefix < 2)
		rc_p->Encode(prefix);
	else if (prefix < 10)
	{
		rc_p->Encode(prefix);
		uint64_t ctx_suf = context_suffix_flag;
		ctx_suf += ((context_t)symbol) << 8;
		ctx_suf += (context_t)prefix;
		uint32_t max_value_for_this_prefix = 1u << (prefix - 1);

		switch (max_value_for_this_prefix)
		{
		case 2:
			find_rce_coder(rce_coders_rl_suf2, ctx_suf)->Encode(len - max_value_for_this_prefix);
			break;
		case 4:
			find_rce_coder(rce_coders_rl_suf4, ctx_suf)->Encode(len - max_value_for_this_prefix);
			break;
		case 8:
			find_rce_coder(rce_coders_rl_suf8, ctx_suf)->Encode(len - max_value_for_this_prefix);
			break;
		case 16:
			find_rce_coder(rce_coders_rl_suf16, ctx_suf)->Encode(len - max_value_for_this_prefix);
			break;
		case 32:
			find_rce_coder(rce_coders_rl_suf32, ctx_suf)->Encode(len - max_value_for_this_prefix);
			break;
		case 64:
			find_rce_coder(rce_coders_rl_suf64, ctx_suf)->Encode(len - max_value_for_this_prefix);
			break;
		case 128:
			find_rce_coder(rce_coders_rl_suf128, ctx_suf)->Encode(len - max_value_for_this_prefix);
			break;
		case 256:
			find_rce_coder(rce_coders_rl_suf256, ctx_suf)->Encode(len - max_value_for_this_prefix);
			break;
		}
	}
	else
	{
		rc_p->Encode(10);		// flag for large value

		context_t ctx_large1 = context_large_value1_flag;
		ctx_large1 += ((context_t)symbol) << 16;
		auto rc_l1 = find_rce_coder(rce_coders_large_val, ctx_large1);
		uint32_t lv1 = (len >> 16) & 0xff;
		rc_l1->Encode(lv1);

		context_t ctx_large2 = context_large_value2_flag;
		ctx_large2 += ((context_t)symbol) << 16;
		ctx_large2 += (context_t)lv1;
		auto rc_l2 = find_rce_coder(rce_coders_large_val, ctx_large2);
		uint32_t lv2 = (len >> 8) & 0xff;
		rc_l2->Encode(lv2);

		context_t ctx_large3 = context_large_value3_flag;
		ctx_large3 += ((context_t)symbol) << 16;
		ctx_large3 += ((context_t)lv1) << 8;
		ctx_large3 += (context_t)lv2;
		auto rc_l3 = find_rce_coder(rce_coders_large_val, ctx_large3);
		uint32_t lv3 = len & 0xff;
		rc_l3->Encode(lv3);
	}
}

// ************************************************************************************
void CCompressedFile::decode_run_len(uint32_t& symbol, uint32_t& len)
{
	// Decode symbol
	auto rc_sym = find_rcd_coder(rcd_coders_rl_sym, ctx_symbol + context_symbol_flag);
	symbol = (uint8_t)rc_sym->Decode();

	if (symbol == 15)
	{
		while (true)
		{
			uint32_t x = (uint8_t)rc_sym->Decode();

			symbol += x;

			if (x < 15)
				break;
		}
	}

	uint32_t symbol_normalized = (symbol > 15) ? 15 : symbol;

	ctx_symbol <<= 4;
	ctx_symbol += (context_t)symbol_normalized;
	ctx_symbol &= context_symbol_mask;

	rcd_coders_rl_sym.prefetch(ctx_symbol + context_symbol_flag);

	ctx_prefix <<= 4;
	ctx_prefix += (context_t)symbol_normalized;
	ctx_prefix &= context_prefix_mask;

	// Decode run length
	auto rc_p = find_rcd_coder(rcd_coders_rl_pref, ctx_prefix + context_prefix_flag);

	uint32_t prefix = rc_p->Decode();

	if (prefix < 2)
		len = prefix;
	else if (prefix < 10)
	{
		uint64_t ctx_suf = context_suffix_flag;
		ctx_suf += ((context_t)symbol_normalized) << 8;
		ctx_suf += (context_t)prefix;
		uint32_t max_value_for_this_prefix = 1u << (prefix - 1);

		switch (max_value_for_this_prefix)
		{
		case 2:
			len = max_value_for_this_prefix + find_rcd_coder(rcd_coders_rl_suf2, ctx_suf)->Decode();
			break;
		case 4:
			len = max_value_for_this_prefix + find_rcd_coder(rcd_coders_rl_suf4, ctx_suf)->Decode();
			break;
		case 8:
			len = max_value_for_this_prefix + find_rcd_coder(rcd_coders_rl_suf8, ctx_suf)->Decode();
			break;
		case 16:
			len = max_value_for_this_prefix + find_rcd_coder(rcd_coders_rl_suf16, ctx_suf)->Decode();
			break;
		case 32:
			len = max_value_for_this_prefix + find_rcd_coder(rcd_coders_rl_suf32, ctx_suf)->Decode();
			break;
		case 64:
			len = max_value_for_this_prefix + find_rcd_coder(rcd_coders_rl_suf64, ctx_suf)->Decode();
			break;
		case 128:
			len = max_value_for_this_prefix + find_rcd_coder(rcd_coders_rl_suf128, ctx_suf)->Decode();
			break;
		case 256:
			len = max_value_for_this_prefix + find_rcd_coder(rcd_coders_rl_suf256, ctx_suf)->Decode();
			break;
		}
	}
	else
	{
		context_t ctx_large1 = context_large_value1_flag;
		ctx_large1 += ((context_t)symbol_normalized) << 16;
		auto rc_l1 = find_rcd_coder(rcd_coders_large_val, ctx_large1);
		uint32_t lv1 = rc_l1->Decode();

		context_t ctx_large2 = context_large_value2_flag;
		ctx_large2 += ((context_t)symbol_normalized) << 16;
		ctx_large2 += (context_t)lv1;
		auto rc_l2 = find_rcd_coder(rcd_coders_large_val, ctx_large2);
		uint32_t lv2 = rc_l2->Decode();

		context_t ctx_large3 = context_large_value3_flag;
		ctx_large3 += ((context_t)symbol_normalized) << 16;
		ctx_large3 += ((context_t)lv1) << 8;
		ctx_large3 += (context_t)lv2;
		auto rc_l3 = find_rcd_coder(rcd_coders_large_val, ctx_large3);
		uint32_t lv3 = rc_l3->Decode();

		len = (lv1 << 16) + (lv2 << 8) + lv3;

		prefix = ilog2(len);
	}

	ctx_prefix <<= 4;
	ctx_prefix += (context_t)prefix;
	ctx_prefix &= context_prefix_mask;

	rcd_coders_rl_pref.prefetch(ctx_prefix + context_prefix_flag);
}

// ******************************************************************************
bool CCompressedFile::OptimizeDB(function_size_graph_t& _function_size_graph, function_data_graph_t& _function_data_graph)
{
	cout << "Archive optimization\n";

	function_size_graph = _function_size_graph;
	function_data_graph = _function_data_graph;

	if (tmp_archive)
		delete tmp_archive;

	tmp_archive = new CArchive(true);

	string tmp_name = archive_name + "_vcfshark_tmp";
	rename(archive_name.c_str(), tmp_name.c_str());

	if (!tmp_archive->Open(tmp_name) || !archive->Open(archive_name))
	{
		std::cerr << "Cannot open archive\n";
		exit(1);
	}

	const vector<string> meta_stream_names = {
		"db_chrom_size", "db_pos_size", "db_id_size", "db_ref_size", "db_alt_size", "db_qual_size",
		"idb_chrom_data", "idb_pos_data", "idb_id_data", "idb_ref_data", "idb_alt_data", "idb_qual_data",
		"db_params", "db_meta", "db_header", "db_samples"
	};

	//int no_keys = (tmp_archive->GetNoStreams() - meta_stream_names.size()) / 2;
	vector<uint8_t> vec;
//	size_t meta;

	process_function_size(no_keys, v_size_nodes, v_size_edges);
//	process_function_data(no_keys, v_data_nodes, v_data_edges);
	process_function_data_eq_only(no_keys, v_data_nodes, v_data_edges);

	// Store description of size and data
	store_nodes("size_nodes", v_size_nodes);
	store_edges("size_edges", v_size_edges, (int) v_size_nodes.size());
	store_nodes("data_nodes", v_data_nodes);
	store_edges("data_edges", v_data_edges, (int) v_data_nodes.size());

	// Process key fields
	for (uint32_t i = 0; i < no_keys; ++i)
		if (v_size_nodes[i].second)
			copy_stream("key_" + to_string(v_size_nodes[i].first) + "_size");
		else
		{
			pair<int, int> pid;

			for (auto& x : v_size_edges)
				if (x.second == v_size_nodes[i].first)
					pid = x;

			link_stream("key_" + to_string(v_size_nodes[i].first) + "_size", "key_" + to_string(pid.first) + "_size");
		}

	for (uint32_t i = 0; i < no_keys; ++i)
		if (v_data_nodes[i].second)
			copy_stream("key_" + to_string(v_data_nodes[i].first) + "_data");
		else
		{
			pair<int, int> pid;

			for (auto& x : v_data_edges)
				if (x.second == v_data_nodes[i].first)
					pid = x;

			link_stream("key_" + to_string(v_data_nodes[i].first) + "_data", "key_" + to_string(pid.first) + "_data");
//			store_function("func_" + to_string(v_data_nodes[i].first) + "_data", pid.first, function_data_graph[pid]);
		}

	for (auto sn : meta_stream_names)
		copy_stream(sn);

	tmp_archive->Close();
	archive->Close();
	remove(tmp_name.c_str());
	cout << endl;

	return true;
}

// ******************************************************************************
bool CCompressedFile::process_function_size(int no_keys,
	vector<pair<int, bool>>& v_out_nodes, vector<pair<int, int>>& v_out_edges)
{
	vector<CGraphOptimizer::node_t> v_in_nodes;
	vector<CGraphOptimizer::edge_t> v_in_edges;

	unordered_map<uint64_t, pair<int, uint64_t>> node_hashes;
	vector<uint8_t> v_data;
	vector<uint8_t> v_data_src;
	size_t metadata, metadata_src;

	for (int i = 0; i < no_keys; ++i)
	{
		string ks = "key_" + to_string(i) + "_size";
		auto iks = tmp_archive->GetStreamId(ks);

		uint64_t h = 0;
		uint64_t s_size = 0;

		while (tmp_archive->GetPart(iks, v_data, metadata))
		{
			for (auto c : v_data)
				h += ((uint64_t)c) * 127ull;
			s_size += v_data.size();
		}

		auto p = node_hashes.find(h);
		
		if (p != node_hashes.end() && p->second.second == s_size)
		{
			auto iks_src = tmp_archive->GetStreamId("key_" + to_string(p->second.first) + "_size");
			bool same = true;

			tmp_archive->ResetStreamPartIterator(iks);
			tmp_archive->ResetStreamPartIterator(iks_src);

			while (tmp_archive->GetPart(iks, v_data, metadata) && tmp_archive->GetPart(iks_src, v_data_src, metadata_src) && same)
				same = v_data == v_data_src && metadata == metadata_src;

			if (same)
			{
				v_out_edges.emplace_back(p->second.first, i);
				v_out_nodes.emplace_back(i, false);
				continue;
			}
		}

		node_hashes[h] = make_pair(i, s_size);
		v_out_nodes.emplace_back(i, true);
	}
	
	return true;
}

// ******************************************************************************
bool CCompressedFile::process_function_data(int no_keys,
	vector<pair<int, bool>>& v_out_nodes, vector<pair<int, int>>& v_out_edges)
{
	vector<CGraphOptimizer::node_t> v_in_nodes;
	vector<CGraphOptimizer::edge_t> v_in_edges;

	for (int i = 0; i < no_keys; ++i)
	{
		string ks = "key_" + to_string(i) + "_data";
		auto iks = tmp_archive->GetStreamId(ks);

		v_in_nodes.emplace_back(i, tmp_archive->GetCompressedSize(iks));
	}

	for (auto& edge : function_data_graph)
	{
		uint64_t cost = 0;
		bool equality = true;

		for (auto e : edge.second)
		{
			cost += e.first.size() + 1 + e.second.size() + 1;
			equality &= e.first == e.second;
		}

		v_in_edges.emplace_back(edge.first.first, edge.first.second, equality, cost);
	}

	CGraphOptimizer graph_optimizer;

	graph_optimizer.Optimize(v_in_nodes, v_in_edges, v_out_nodes, v_out_edges);

	return true;
}

// ******************************************************************************
bool CCompressedFile::process_function_data_eq_only(int no_keys,
	vector<pair<int, bool>>& v_out_nodes, vector<pair<int, int>>& v_out_edges)
{
	vector<CGraphOptimizer::node_t> v_in_nodes;
	vector<CGraphOptimizer::edge_t> v_in_edges;

	unordered_map<uint64_t, pair<int, uint64_t>> node_hashes;
	vector<uint8_t> v_data;
	vector<uint8_t> v_data_src;
	size_t metadata, metadata_src;

	for (int i = 0; i < no_keys; ++i)
	{
		string ks = "key_" + to_string(i) + "_data";
		auto iks = tmp_archive->GetStreamId(ks);

		uint64_t h = 0;
		uint64_t s_size = 0;

		while (tmp_archive->GetPart(iks, v_data, metadata))
		{
			for (auto c : v_data)
				h += ((uint64_t)c) * 127ull;
			s_size += v_data.size();
		}

		auto p = node_hashes.find(h);

		if (p != node_hashes.end() && p->second.second == s_size)
		{
			auto iks_src = tmp_archive->GetStreamId("key_" + to_string(p->second.first) + "_data");
			bool same = true;

			tmp_archive->ResetStreamPartIterator(iks);
			tmp_archive->ResetStreamPartIterator(iks_src);

			while (tmp_archive->GetPart(iks, v_data, metadata) && tmp_archive->GetPart(iks_src, v_data_src, metadata_src) && same)
				same = v_data == v_data_src && metadata == metadata_src;

			if (same)
			{
				v_out_edges.emplace_back(p->second.first, i);
				v_out_nodes.emplace_back(i, false);
				continue;
			}
		}

		node_hashes[h] = make_pair(i, s_size);
		v_out_nodes.emplace_back(i, true);
	}

	return true;
}

// ******************************************************************************
void CCompressedFile::store_nodes(string stream_name, vector<pair<int, bool>>& v_nodes)
{
	auto sid = archive->RegisterStream(stream_name);

	vector<uint8_t> vec;
	uint32_t nb = (uint32_t) no_bytes(no_keys);

	for (auto x : v_nodes)
	{
		vec.emplace_back((uint8_t)x.second);
		for (uint32_t i = 0; i < nb; ++i)
			vec.emplace_back((x.first >> (8 * i)) & 0xff);
	}

	archive->AddPart(sid, vec, v_nodes.size());
}

// ******************************************************************************
void CCompressedFile::load_nodes(string stream_name, vector<pair<int, bool>>& v_nodes)
{
	auto sid = archive->GetStreamId(stream_name);

	v_nodes.clear();

	vector<uint8_t> vec;
	uint32_t nb = (uint32_t) no_bytes(no_keys);
	size_t metadata;

	v_nodes.resize(no_keys);

	archive->GetPart(sid, vec, metadata);
	auto p = vec.begin();

	for (auto& x : v_nodes)
	{
		x.second = (bool)*p++;

		x.first = 0;
		for (uint32_t i = 0; i < nb; ++i)
			x.first += ((int)*p++) << (8 * i);
	}
}

// ******************************************************************************
void CCompressedFile::store_edges(string stream_name, vector<pair<int, int>>& v_edges, int no_keys)
{
	auto sid = archive->RegisterStream(stream_name);

	vector<uint8_t> vec;
	uint32_t nb = (uint32_t) no_bytes(no_keys);

	for (auto x : v_edges)
	{
		for (uint32_t i = 0; i < nb; ++i)
			vec.emplace_back((x.first >> (8 * i)) & 0xff);
		for (uint32_t i = 0; i < nb; ++i)
			vec.emplace_back((x.second >> (8 * i)) & 0xff);
	}

	archive->AddPart(sid, vec, v_edges.size());
}

// ******************************************************************************
void CCompressedFile::load_edges(string stream_name, vector<pair<int, int>>& v_edges, int no_keys)
{
	auto sid = archive->GetStreamId(stream_name);

	vector<uint8_t> vec;
	uint32_t nb = (uint32_t) no_bytes(no_keys);
	size_t metadata;

	archive->GetPart(sid, vec, metadata);
	auto p = vec.begin();

	v_edges.resize(metadata);

	for (auto& x : v_edges)
	{
		x.first = 0;
		for (uint32_t i = 0; i < nb; ++i)
			x.first += ((int)*p++) << (8 * i);

		x.second = 0;
		for (uint32_t i = 0; i < nb; ++i)
			x.second += ((int)*p++) << (8 * i);
	}
}

// ******************************************************************************
void CCompressedFile::copy_stream(string stream_name)
{
	auto in_iks = tmp_archive->GetStreamId(stream_name);
	auto out_iks = archive->RegisterStream(stream_name);
	vector<uint8_t> vec;
	size_t meta;

	tmp_archive->ResetStreamPartIterator(in_iks);
	while (tmp_archive->GetPart(in_iks, vec, meta))
		archive->AddPart(out_iks, vec, meta);

	archive->SetRawSize(out_iks, tmp_archive->GetRawSize(in_iks));
}

// ******************************************************************************
void CCompressedFile::link_stream(string stream_name, string target_name)
{
	auto out_iks = archive->RegisterStream(stream_name);
	auto target_iks = archive->GetStreamId(target_name);

	archive->LinkStream(out_iks, stream_name, target_iks);
}

// ******************************************************************************
void CCompressedFile::store_function(string stream_name, int src_id, function_data_item_t& func)
{
	bool equality = true;
	vector<uint8_t> vec;

	append(vec, src_id);

	for (auto& x : func)
		if (x.first != x.second)
		{
			equality = false;
			break;
		}

	if (equality)
		vec.emplace_back(1);
	else
	{
		vec.emplace_back(0);

		for (auto& x : func)
		{
			vec.emplace_back((uint8_t)x.first.size());
			for (auto c : x.first)
				vec.emplace_back(c);

			vec.emplace_back((uint8_t)x.second.size());
			for (auto c : x.second)
				vec.emplace_back(c);
		}
	}

	auto iks = archive->RegisterStream(stream_name);
	archive->AddPart(iks, vec, func.size());
	archive->SetRawSize(iks, 0);
}

// ******************************************************************************
void CCompressedFile::load_function(string stream_name, int &src_id, function_data_item_t& func)
{
	vector<uint8_t> vec;
	size_t metadata;

	auto iks = archive->GetStreamId(stream_name);
	archive->GetPart(iks, vec, metadata);

	if (vec.empty())
		return;			// equality

	auto p = vec.begin();
	func.clear();
	
	size_t offset = 0;
	int64_t src_id64;

	read(vec, offset, src_id64);
	p += offset;

	if (*p++ == 1)		// equality 
		return;
	else
	{
		vector<uint8_t> v_from, v_to;

		for (uint32_t i = 0; i < metadata; ++i)
		{
			v_from.resize(*p++);
			for (auto& x : v_from)
				x = *p++;

			v_to.resize(*p++);
			for (auto& x : v_to)
				x = *p++;

			func[v_from] = v_to;
		}
	}
}

// ******************************************************************************
void CCompressedFile::store_function(string stream_name, int src_id, function_size_item_t& func)
{
	bool equality = true;
	vector<uint8_t> vec;

	append(vec, src_id);

	for (auto& x : func)
		if (x.first != x.second)
		{
			equality = false;
			break;
		}

	if (equality)
		vec.emplace_back(1);
	else
	{
		vec.emplace_back(0);

		for (auto& x : func)
		{
			auto nb1 = no_bytes(x.first);
			auto nb2 = no_bytes(x.second);

			vec.emplace_back(nb1);
			for (uint32_t i = 0; i < nb1; ++i)
				vec.emplace_back((x.first >> (8 * i)) & 0xff);

			vec.emplace_back(nb2);
			for (uint32_t i = 0; i < nb2; ++i)
				vec.emplace_back((x.second >> (8 * i)) & 0xff);
		}
	}

	auto iks = archive->RegisterStream(stream_name);
	archive->AddPart(iks, vec, func.size());
	archive->SetRawSize(iks, 0);
}

// ******************************************************************************
void CCompressedFile::load_function(string stream_name, int &src_id, function_size_item_t& func)
{
	vector<uint8_t> vec;
	size_t metadata;

	auto iks = archive->GetStreamId(stream_name);
	archive->GetPart(iks, vec, metadata);

	auto p = vec.begin();
	func.clear();

	size_t offset = 0;
	int64_t src_id64;

	read(vec, offset, src_id64);
	p += offset;

	if (*p++ == 1)		// equality 
		return;
	else
	{
		uint32_t v_from, v_to;

		for (uint32_t i = 0; i < metadata; ++i)
		{
			v_from = v_to = 0;

			uint32_t nb1 = *p++;
			for (uint32_t i = 0; i < nb1; ++i)
				v_from += ((uint32_t)*p++) << (8 * i);

			uint32_t nb2 = *p++;
			for (uint32_t i = 0; i < nb2; ++i)
				v_to += ((uint32_t)*p++) << (8 * i);

			func[v_from] = v_to;
		}
	}
}

// EOF
