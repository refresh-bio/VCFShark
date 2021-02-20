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
#include <functional>
#include <vector>

using namespace std;

#include "cfile.h"
#include "utils.h"

// ************************************************************************************
CCompressedFile::CCompressedFile()
{
	open_mode = open_mode_t::none;

	vcs_compression_level = 3;

	archive = nullptr;
	tmp_archive = nullptr;

	q_packages = nullptr;
	q_preparation_ids = nullptr;

	vios_i = new CVectorIOStream(v_vios_i);
	vios_o = new CVectorIOStream(v_vios_o);

	rce = nullptr;
	rcd = nullptr;
}

// ************************************************************************************
CCompressedFile::~CCompressedFile()
{
	v_coder_threads.clear();

	using fo_t = function<void(void)>;

	CRegisteringQueue<fo_t> q_fo(1);

	vector<thread> v_thr;

	v_thr.reserve(no_coder_threads);

	for (uint32_t i = 0; i < no_coder_threads; ++i)
		v_thr.emplace_back(thread([&, this]() {

		fo_t fo;

		while (!q_fo.IsCompleted())
		{
			if (q_fo.Pop(fo))
				fo();
		}
	}));

	if (rce)
		q_fo.Push([=] {delete rce; });
	if (rcd)
		q_fo.Push([=] {delete rcd; });

	if (vios_i)
		q_fo.Push([=] {delete vios_i; });
	if (vios_o)
		q_fo.Push([=] {delete vios_o; });

	if (q_packages)
		q_fo.Push([=] {delete q_packages; });

	if (q_preparation_ids)
		q_fo.Push([=] {delete q_preparation_ids; });

	for (auto p : v_bsc_size)
		if (p)
			q_fo.Push([=] {delete p; });

	for (auto p : v_bsc_data)
		if (p)
			q_fo.Push([=] {	delete p; });

	for (auto p : v_bsc_db_size)
		if (p)
			q_fo.Push([=] {delete p; });

	for (auto p : v_bsc_db_data)
		if (p)
			q_fo.Push([=] {delete p; });

	for (auto p : v_format_compress)
		if(p)
			q_fo.Push([=] {delete p; });

	for (auto p : v_packages)
		if (p)
			q_fo.Push([=] {delete p;});

	for (auto p : v_db_packages)
		if (p)
			q_fo.Push([=] {delete p;});

	if (archive)
		q_fo.Push([=] {delete archive; });

	if (tmp_archive)
		q_fo.Push([=] {delete tmp_archive; });
		
	q_fo.MarkCompleted();

	for (auto& t : v_thr)
		t.join();
}

// ************************************************************************************
bool CCompressedFile::OpenForReading(string file_name)
{
	prev_pos = 0;
	archive_name = file_name;
	archive = new CArchive(true);

	CBSCWrapper::InitLibrary(p_bsc_features);

	if (!archive->Open(file_name))
	{
		cerr << "Cannot open " << file_name << "\n";
		return false;
	}

    load_descriptions();

	load_nodes("size_nodes", v_size_nodes);
	load_edges("size_edges", v_size_edges, (int) v_size_nodes.size());
	load_nodes("data_nodes", v_data_nodes);
	load_edges("data_edges", v_data_edges, (int) v_data_nodes.size());

	gt_stream_id = archive->GetStreamId("key_" + to_string(gt_key_id) + "_size");

	q_preparation_ids = new CRegisteringQueue<pair<int, int>>(1);

	m_data_nodes.clear();
	m_data_nodes.resize(no_keys, true);

	m_data_edges.clear();
	m_data_edges.resize(no_keys);
	for(auto e : v_data_edges)
		m_data_edges[e.second] = e.first;

	v_packages.resize(no_keys, nullptr);
	for (uint32_t i = 0; i < no_keys; ++i)
		q_preparation_ids->Push(make_pair(i, -1));

	v_db_packages.resize(no_db_fields, nullptr);
	for(uint32_t i = 0; i < no_db_fields; ++i)
		q_preparation_ids->Push(make_pair(-1, i));

	v_i_buf.resize(no_keys);
	v_i_db_buf.resize(no_db_fields);

	open_mode = open_mode_t::reading;

	v_coder_threads.reserve(no_coder_threads);

	v_bsc_size.resize(no_keys);
	v_bsc_data.resize(no_keys);
	v_text_pp.resize(no_keys);
	v_bsc_db_size.resize(no_db_fields);
	v_bsc_db_data.resize(no_db_fields);

	v_format_compress.resize(no_keys, nullptr);

	rcd = new CRangeDecoder<CVectorIOStream>(*vios_i);

	pbwt_initialised = false;

	InitPBWT();

	for (uint32_t i = 0; i < no_keys; ++i)
	{
		v_bsc_size[i] = new CBSCWrapper;
		v_bsc_size[i]->InitDecompress();

		v_bsc_data[i] = new CBSCWrapper;

		if (keys[i].keys_type == key_type_t::fmt || keys[i].keys_type == key_type_t::info)
		{
			v_format_compress[i] = new CFormatCompress("key " + to_string(i), vcs_compression_level);
			v_format_compress[i]->SetNoSamples(no_samples);
		}

		switch (keys[i].type)
		{
		case BCF_HT_FLAG:
			v_bsc_data[i]->InitDecompress();
			break;
		case BCF_HT_INT:
			v_bsc_data[i]->InitDecompress();
			break;
		case BCF_HT_REAL:
			v_bsc_data[i]->InitDecompress();
			break;
		case BCF_HT_STR:
			v_bsc_data[i]->InitDecompress();
			break;
		}
	}

	for (uint32_t i = 0; i < no_db_fields; ++i)
	{
		v_bsc_db_size[i] = new CBSCWrapper;
		v_bsc_db_size[i]->InitDecompress();
		v_bsc_db_data[i] = new CBSCWrapper;
		v_bsc_db_data[i]->InitDecompress();
	}

	for (uint32_t i = 0; i < no_coder_threads; ++i)
		v_coder_threads.emplace_back(thread([&]() {

		while (!q_preparation_ids->IsCompleted())
		{
			SPackage* pck = new SPackage;
			size_t raw_size;
			vector<uint8_t> v_tmp;
			pair<int, int> p_ids;

			if (!q_preparation_ids->Pop(p_ids))
			{
				delete pck;
				break;
			}

			if (p_ids.first >= 0)
			{
				pck->stream_id_size = archive->GetStreamId("key_" + to_string(p_ids.first) + "_size");

                pck->is_func = !m_data_nodes[p_ids.first];

				if (archive->GetPart(pck->stream_id_size, pck->v_compressed, raw_size))
				{
					if ((int) pck->stream_id_size != gt_stream_id)		// keys
					{
						pck->key_id = p_ids.first;
						if (keys[pck->key_id].keys_type == key_type_t::fmt && keys[pck->key_id].type != BCF_HT_STR)
							decompress_format(pck, raw_size, v_tmp);
						else if (keys[pck->key_id].keys_type == key_type_t::info &&
							(keys[pck->key_id].type == BCF_HT_INT || keys[pck->key_id].type == BCF_HT_REAL))
							decompress_info(pck, raw_size, v_tmp);
						else
							decompress_field(pck, raw_size, v_tmp);
					}
					else
					{
						pck->key_id = p_ids.first;
						decompress_gt(pck, raw_size);
					}

					lock_guard<mutex> lck(m_packages);
					v_packages[pck->key_id] = pck;
				}
				else
				{
					pck->v_size.clear();
					pck->v_data.clear();

					lock_guard<mutex> lck(m_packages);

					pck->key_id = p_ids.first;
					v_packages[pck->key_id] = pck;
				}
			}
			else
			{
				pck->db_id = p_ids.second;
				pck->stream_id_size = archive->GetStreamId(db_stream_name_size[p_ids.second]);
				pck->stream_id_data = archive->GetStreamId(db_stream_name_data[p_ids.second]);
				
				if (archive->GetPart(pck->stream_id_size, pck->v_compressed, raw_size))
				{
					decompress_db(pck, raw_size, v_tmp);
					lock_guard<mutex> lck(m_packages);
					v_db_packages[pck->db_id] = pck;
				}
				else
				{
					pck->v_size.clear();
					pck->v_data.clear();

					lock_guard<mutex> lck(m_packages);
					v_db_packages[pck->db_id] = pck;
				}
			}
						
			cv_packages.notify_all();
		}
			}));

	return true;
}

// ************************************************************************************
bool CCompressedFile::OpenForWriting(string file_name, uint32_t _no_keys)
{
	prev_pos = 0;
	archive_name = file_name;
	if (archive)
		delete archive;
	archive = new CArchive(false);

	CBSCWrapper::InitLibrary(p_bsc_features);

	if (!archive->Open(file_name))
	{
		cerr << "Cannot open " << file_name << "\n";
		return false;
	}

    no_keys = _no_keys;
    
	v_o_buf.resize(no_keys);
	v_o_db_buf.resize(no_db_fields);

	v_buf_ids_size.resize(no_keys, -1);
	v_buf_ids_data.resize(no_keys, -1);
	v_buf_ids_func.resize(no_keys, -1);

	v_bsc_size.resize(no_keys);
	v_bsc_data.resize(no_keys);
	v_text_pp.resize(no_keys);
	v_coder_part_ids.resize(no_keys + no_db_fields, 0);
	v_text_part_ids.resize(no_keys + no_db_fields, 0);

	v_format_compress.resize(no_keys, nullptr);

	open_mode = open_mode_t::writing;
	pbwt_initialised = false;
	no_variants = 0;

	InitPBWT();

	for (uint32_t i = 0; i < no_keys; i++)
	{
		if((int) i != gt_key_id)
			v_o_buf[i].SetMaxSize(max_buffer_size, i * max_buffer_size / no_keys);
		else
			v_o_buf[i].SetMaxSize(max_buffer_gt_size, 0);

		v_bsc_size[i] = new CBSCWrapper;
		v_bsc_size[i]->InitCompress(p_bsc_size);
		v_bsc_data[i] = new CBSCWrapper;

		v_buf_ids_size[i] = archive->RegisterStream("key_" + to_string(i) + "_size");

		if (keys[i].keys_type == key_type_t::fmt || keys[i].keys_type == key_type_t::info)
		{
			v_format_compress[i] = new CFormatCompress("key " + to_string(i), vcs_compression_level);
			v_format_compress[i]->SetNoSamples(no_samples);
		}

		switch (keys[i].type)
		{
		case BCF_HT_FLAG:
			v_bsc_data[i]->InitCompress(p_bsc_flag);
			break;
		case BCF_HT_INT:
			v_bsc_data[i]->InitCompress(p_bsc_int);
			break;
		case BCF_HT_REAL:
			v_bsc_data[i]->InitCompress(p_bsc_real);
			break;
		case BCF_HT_STR:
			v_bsc_data[i]->InitCompress(p_bsc_text);
			break;
		}
	}

	if (rce)
		delete rce;
	rce = new CRangeEncoder<CVectorIOStream>(*vios_o);

	for (uint32_t i = 0; i < no_keys; i++)
		v_buf_ids_data[i] = archive->RegisterStream("key_" + to_string(i) + "_data");

	// Register streams for variant descriptions
	v_db_ids_size.clear();
	for(auto x : db_stream_name_size)
		v_db_ids_size.emplace_back(archive->RegisterStream(x));

	v_db_ids_data.clear();
	for (auto x : db_stream_name_data)
		v_db_ids_data.emplace_back(archive->RegisterStream(x));

	for(uint32_t i = 0; i < no_db_fields; ++i)
		v_o_db_buf[i].SetMaxSize(max_buffer_db_size, i * max_buffer_db_size / no_db_fields / 8);

	v_bsc_db_size.resize(no_db_fields);
	v_bsc_db_data.resize(no_db_fields);

	v_bsc_db_size[id_db_chrom] = new CBSCWrapper;
	v_bsc_db_size[id_db_chrom]->InitCompress(p_bsc_size);
	v_bsc_db_data[id_db_chrom] = new CBSCWrapper;
	v_bsc_db_data[id_db_chrom]->InitCompress(p_bsc_data);

	v_bsc_db_size[id_db_pos] = new CBSCWrapper;
	v_bsc_db_size[id_db_pos]->InitCompress(p_bsc_size);
	v_bsc_db_data[id_db_pos] = new CBSCWrapper;
	v_bsc_db_data[id_db_pos]->InitCompress(p_bsc_db_pos);

	v_bsc_db_size[id_db_id] = new CBSCWrapper;
	v_bsc_db_size[id_db_id]->InitCompress(p_bsc_size);
	v_bsc_db_data[id_db_id] = new CBSCWrapper;
	v_bsc_db_data[id_db_id]->InitCompress(p_bsc_db_id);

	v_bsc_db_size[id_db_ref] = new CBSCWrapper;
	v_bsc_db_size[id_db_ref]->InitCompress(p_bsc_size);
	v_bsc_db_data[id_db_ref] = new CBSCWrapper;
	v_bsc_db_data[id_db_ref]->InitCompress(p_bsc_db_ref);

	v_bsc_db_size[id_db_alt] = new CBSCWrapper;
	v_bsc_db_size[id_db_alt]->InitCompress(p_bsc_size);
	v_bsc_db_data[id_db_alt] = new CBSCWrapper;
	v_bsc_db_data[id_db_alt]->InitCompress(p_bsc_db_alt);

	v_bsc_db_size[id_db_qual] = new CBSCWrapper;
	v_bsc_db_size[id_db_qual]->InitCompress(p_bsc_size);
	v_bsc_db_data[id_db_qual] = new CBSCWrapper;
	v_bsc_db_data[id_db_qual]->InitCompress(p_bsc_db_qual);

	if (q_packages)
		delete q_packages;
	q_packages = new CRegisteringQueue<SPackage>(1);

	v_cnt_packages.resize(no_keys, 0);
	v_cnt_db_packages.resize(no_db_fields, 0);

	v_coder_threads.reserve(no_coder_threads);

	for (uint32_t i = 0; i < no_coder_threads; ++i)
		v_coder_threads.emplace_back(thread([&]() {
		
		SPackage pck;
		vector<uint8_t> v_compressed;
		vector<uint8_t> v_tmp;

		auto fo = [this](SPackage& pck)->bool {return check_coder_compressor(pck); };

		while (!q_packages->IsCompleted())
		{
			if (!q_packages->PopWithHint<SPackage>(pck, fo))
				continue;

			{
				unique_lock<mutex> lck(m_packages);

				if(pck.type == SPackage::package_t::db)
					--v_cnt_db_packages.at(pck.db_id);
				else
					--v_cnt_packages.at(pck.key_id);
				cv_packages.notify_all();
			}

			if (pck.type == SPackage::package_t::fields)
			{
				if(keys[pck.key_id].keys_type == key_type_t::fmt && keys[pck.key_id].type != BCF_HT_STR)
					compress_format(pck, v_compressed, v_tmp);
//				else if (keys[pck.stream_id].keys_type == key_type_t::info && keys[pck.stream_id].type != BCF_HT_STR)
				else if (keys[pck.key_id].keys_type == key_type_t::info && 
					(keys[pck.key_id].type == BCF_HT_INT || keys[pck.key_id].type == BCF_HT_REAL))
				{
//					compress_format(pck, v_compressed, v_tmp);
					compress_info(pck, v_compressed, v_tmp);
				}
				else
					compress_field(pck, v_compressed, v_tmp);
			}
			else if (pck.type == SPackage::package_t::gt)
				compress_gt(pck);
			else
				compress_db(pck, v_compressed, v_tmp);
		}
			}));

	return true;
}

// ************************************************************************************
bool CCompressedFile::Close()
{
#ifdef LOG_INFO
	cout << "Distinct values:\n";

	for (int i = 0; i < no_keys; ++i)
		cout << i << ": " << distinct_values[i].size() << endl;
#endif

	if (open_mode == open_mode_t::writing)
	{
		vector<uint32_t> v_size;
		vector<uint8_t> v_data;
		vector<uint8_t> v_aux;

		for (uint32_t i = 0; i < no_keys; ++i)
		{
			v_o_buf[i].GetBuffer(v_size, v_data);
			auto part_id = archive->AddPartPrepare(v_buf_ids_size[i]);
			archive->AddPartPrepare(v_buf_ids_data[i]);

			SPackage pck((int) i != gt_key_id ? SPackage::package_t::fields : SPackage::package_t::gt, i, -1, v_buf_ids_size[i], v_buf_ids_data[i], part_id, v_size, v_data, v_aux);

			q_packages->Push(pck);
		}

		for(uint32_t i = 0; i < no_db_fields; ++i)
		{
			v_o_db_buf[i].GetBuffer(v_size, v_data);
			auto part_id = archive->AddPartPrepare(v_db_ids_size[i]);
			archive->AddPartPrepare(v_db_ids_data[i]);

			SPackage pck(SPackage::package_t::db, -1, i, v_db_ids_size[i], v_db_ids_data[i], part_id, v_size, v_data, v_aux);

			q_packages->Push(pck);
		}

		q_packages->MarkCompleted();

		for (uint32_t i = 0; i < no_coder_threads; ++i)
			v_coder_threads[i].join();

		save_descriptions();

		delete rce;
		rce = nullptr;

		archive->Close();
	}
	else if (open_mode == open_mode_t::reading)
	{
		q_preparation_ids->MarkCompleted();

		for (uint32_t i = 0; i < no_coder_threads; ++i)
			v_coder_threads[i].join();
			
		archive->Close();
	}

	open_mode = open_mode_t::none;
	
	return true;
}

// ************************************************************************************
int CCompressedFile::GetNoSamples()
{
	return no_samples;
}

// ************************************************************************************
void CCompressedFile::SetNoSamples(uint32_t _no_samples)
{
	no_samples = _no_samples;
}

// ************************************************************************************
int CCompressedFile::GetNoVariants()
{
	return no_variants;
}

// ************************************************************************************
bool CCompressedFile::GetMeta(string &_v_meta)
{
	_v_meta = v_meta;

	return true;
}

// ************************************************************************************
bool CCompressedFile::SetMeta(string &_v_meta)
{
	v_meta = _v_meta;

	return true;
}

// ************************************************************************************
bool CCompressedFile::GetHeader(string &_v_header)
{
	_v_header = v_header;

	return true;
}

// ************************************************************************************
bool CCompressedFile::SetHeader(string &_v_header)
{
	v_header = _v_header;

	return true;
}

// ************************************************************************************
bool CCompressedFile::AddSamples(vector<string> &_v_samples)
{
	v_samples = _v_samples;

	return true;
}

// ************************************************************************************
bool CCompressedFile::GetSamples(vector<string> &_v_samples)
{
	_v_samples = v_samples;

	return true;
}
// ************************************************************************************
bool CCompressedFile::SetKeys(vector<key_desc> &_keys)
{
    keys = _keys;
    
    return true;
}
// ************************************************************************************
bool CCompressedFile::GetKeys(vector<key_desc> &_keys)
{
    _keys = keys;
    
    return true;
}
// ************************************************************************************
int CCompressedFile::GetNoKeys()
{
    return no_keys;
}

// ************************************************************************************
void CCompressedFile::SetNoKeys(uint32_t _no_keys)
{
    no_keys = _no_keys;
}

// ************************************************************************************
int CCompressedFile::GetGTId()
{
    return gt_key_id;
}

// ************************************************************************************
void CCompressedFile::SetGTId(uint32_t _gt_key_id)
{
    gt_key_id = _gt_key_id;
}

// ************************************************************************************
int CCompressedFile::GetPloidy()
{
	return ploidy;
}

// ************************************************************************************
void CCompressedFile::SetPloidy(int _ploidy)
{
	ploidy = (uint8_t) _ploidy;
}

// ************************************************************************************
void CCompressedFile::SetNoThreads(int _no_threads)
{
	if (_no_threads > 1)
		no_coder_threads = _no_threads - 1;
	else
		no_coder_threads = 1;
}

// ************************************************************************************
void CCompressedFile::SetCompressionLevel(int _compression_level)
{
	vcs_compression_level = _compression_level;
}

// ************************************************************************************
int CCompressedFile::GetNeglectLimit()
{
	return neglect_limit;
}

// ************************************************************************************
void CCompressedFile::SetNeglectLimit(uint32_t _neglect_limit)
{
	neglect_limit = _neglect_limit;
}

// ************************************************************************************
bool CCompressedFile::Eof()
{
	if (open_mode == open_mode_t::reading)
		return i_variant >= no_variants;

	return false;
}

// ************************************************************************************
bool CCompressedFile::GetVariant(variant_desc_t &desc, vector<field_desc> &fields)
{
	desc.chrom.clear();

	if (i_variant >= no_variants)
		return false;

	int64_t pos;

	for (uint32_t i = 0; i < no_db_fields; ++i)
	{
		if (v_i_db_buf[i].IsEmpty())
		{
			unique_lock<mutex> lck(m_packages);

			cv_packages.wait(lck, [&, this] {return v_db_packages[i] != nullptr; });

			v_i_db_buf[i].SetBuffer(v_db_packages[i]->v_size, v_db_packages[i]->v_data);
			delete v_db_packages[i];
			v_db_packages[i] = nullptr;

			q_preparation_ids->Push(make_pair(-1, i));
		}
	}

	char* str = nullptr;
	uint32_t len = 0;
	
	v_i_db_buf[id_db_chrom].ReadText(str, len);
	desc.chrom = string(str, str + len);
	delete[] str;

	v_i_db_buf[id_db_id].ReadText(str, len);
	desc.id = string(str, str + len);
	delete[] str;

	v_i_db_buf[id_db_ref].ReadText(str, len);
	desc.ref = string(str, str + len);
	delete[] str;

	v_i_db_buf[id_db_alt].ReadText(str, len);
	desc.alt = string(str, str + len);
	delete[] str;

	v_i_db_buf[id_db_qual].ReadText(str, len);
	desc.qual = string(str, str + len);
	delete[] str;

	v_i_db_buf[id_db_pos].ReadInt64(pos);
	pos += prev_pos;
	prev_pos = pos;
	desc.pos = pos;

    // Load and set fields
    for(uint32_t i = 0; i < no_keys; i++)
    {
		int ii = v_data_nodes[i].first;		// Change of column ordering

		if (v_i_buf[ii].IsEmpty())
		{
			unique_lock<mutex> lck(m_packages);

			cv_packages.wait(lck, [&, this] {return v_packages[ii] != nullptr; });

			if (v_packages[ii]->is_func)
				v_i_buf[ii].SetFunction(v_packages[ii]->fun);
			else
				v_i_buf[ii].SetBuffer(v_packages[ii]->v_size, v_packages[ii]->v_data);
			delete v_packages[ii];
			v_packages[ii] = nullptr;

			q_preparation_ids->Push(make_pair(ii, -1));
		}

		switch (keys[ii].type)
		{
		case BCF_HT_INT:
			if(m_data_nodes[ii])
				v_i_buf[ii].ReadInt(fields[ii].data, fields[ii].data_size);
			else
				v_i_buf[ii].FuncInt(fields[ii].data, fields[ii].data_size, fields[m_data_edges[ii]].data, fields[m_data_edges[ii]].data_size);
			fields[ii].present = fields[ii].data != nullptr;
			break;
		case BCF_HT_REAL:
			if (m_data_nodes[ii])
				v_i_buf[ii].ReadReal(fields[ii].data, fields[ii].data_size);
			else
				v_i_buf[ii].FuncReal(fields[ii].data, fields[ii].data_size, fields[m_data_edges[ii]].data, fields[m_data_edges[ii]].data_size);

			fields[ii].present = fields[ii].data != nullptr;
			break;
		case BCF_HT_STR:
			v_i_buf[ii].ReadText(fields[ii].data, fields[ii].data_size);
			fields[ii].present = fields[ii].data != nullptr;
			break;
		case BCF_HT_FLAG:
			uint8_t tmp;
			v_i_buf[ii].ReadFlag(tmp);
			fields[ii].present = (bool)tmp;
			fields[ii].data_size = 1;
			break;
		}
    }
    
	++i_variant;

	return true;
}

// ************************************************************************************
bool CCompressedFile::SetVariant(variant_desc_t &desc, vector<field_desc> &fields)
{
    // Store variant description
	v_o_db_buf[id_db_chrom].WriteText((char*) desc.chrom.c_str(), (uint32_t) desc.chrom.size());
	v_o_db_buf[id_db_pos].WriteInt64(desc.pos - prev_pos);
	v_o_db_buf[id_db_id].WriteText((char*) desc.id.c_str(), (uint32_t) desc.id.size());
	v_o_db_buf[id_db_ref].WriteText((char*) desc.ref.c_str(), (uint32_t) desc.ref.size());
	v_o_db_buf[id_db_alt].WriteText((char*) desc.alt.c_str(), (uint32_t) desc.alt.size());
	v_o_db_buf[id_db_qual].WriteText((char*) desc.qual.c_str(), (uint32_t) desc.qual.size());

	for(uint32_t i = 0; i < no_db_fields; ++i)
		if (v_o_db_buf[i].IsFull())
		{
			auto part_id = archive->AddPartPrepare(v_db_ids_size[i]);
			archive->AddPartPrepare(v_db_ids_data[i]);

			vector<uint32_t> v_size;
			vector<uint8_t> v_data;
			vector<uint8_t> v_aux;

			v_o_db_buf[i].GetBuffer(v_size, v_data);

			SPackage pck(SPackage::package_t::db, -1, i, v_db_ids_size[i], v_db_ids_data[i], part_id, v_size, v_data, v_aux);

			{
				unique_lock<mutex> lck(m_packages);

				cv_packages.wait(lck, [&, this] {return v_cnt_db_packages[i] < max_cnt_packages; });
				++v_cnt_db_packages[i];
			}

			q_packages->Emplace(pck);
		}

	prev_pos = desc.pos;

    for(uint32_t i = 0; i < no_keys; i++)
    {
		switch (keys[i].type)
		{
		case BCF_HT_INT:
			v_o_buf[i].WriteInt(fields[i].data, fields[i].present ? fields[i].data_size : 0);

#ifdef LOG_INFO
			{
				uint32_t* pp = (uint32_t*)fields[i].data;
				for (int j = 0; j < fields[i].data_size; ++j)
					distinct_values[i].insert(pp[j]);
			}
#endif

			break;
		case BCF_HT_REAL:
			v_o_buf[i].WriteReal(fields[i].data, fields[i].present ? fields[i].data_size : 0);

#ifdef LOG_INFO
			{
				uint32_t* pp = (uint32_t*)fields[i].data;

				for (int j = 0; j < fields[i].data_size; ++j)
					distinct_values[i].insert(pp[j]);
			}
#endif

			break;
		case BCF_HT_STR:
			v_o_buf[i].WriteText(fields[i].data, fields[i].present ? fields[i].data_size : 0);
			break;
		case BCF_HT_FLAG:
			v_o_buf[i].WriteFlag(fields[i].present);
			break;
		}

		if (v_o_buf[i].IsFull())
		{
			auto part_id = archive->AddPartPrepare(v_buf_ids_size[i]);
			archive->AddPartPrepare(v_buf_ids_data[i]);

			vector<uint32_t> v_size;
			vector<uint8_t> v_data;
			vector<uint8_t> v_aux;
			
			v_o_buf[i].GetBuffer(v_size, v_data);

			SPackage pck((int) i != gt_key_id ? SPackage::package_t::fields : SPackage::package_t::gt, i, -1, v_buf_ids_size[i], v_buf_ids_data[i], part_id, v_size, v_data, v_aux);

			{
				unique_lock<mutex> lck(m_packages);

				cv_packages.wait(lck, [&, this] {return v_cnt_packages[i] < max_cnt_packages; });
				++v_cnt_packages[i];
			}

			q_packages->Emplace(pck);
		}
    }

	++no_variants;

	return true;
}

// ************************************************************************************
bool CCompressedFile::InitPBWT()
{
	if (open_mode == open_mode_t::reading)
	{
		pbwt.StartReverse(no_samples * ploidy, neglect_limit);
	}
	else if(open_mode == open_mode_t::writing)
	{
		pbwt.StartForward(no_samples * ploidy, neglect_limit);
		pbwt_initialised = true;
	}

	return pbwt_initialised;
}

// EOF
