// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include "application.h"
#include "utils.h"
#include "graph_opt.h"

#include <iostream>
#include <vector>

#include <chrono>
using namespace std::chrono;

using namespace std;

// ******************************************************************************
CApplication::CApplication(const CParams &_params)
{
	params = _params;

	// To avoid unnecessary reallocation of buffers
	v_vcf_data_compress.reserve(no_variants_in_buf);
	v_vcf_data_io.reserve(no_variants_in_buf);
}

// ******************************************************************************
CApplication::~CApplication()
{
}

// ******************************************************************************
bool CApplication::CompressDB()
{
	CBarrier barrier(4);
	unique_ptr<CVCF> vcf(new CVCF());
	unique_ptr<CCompressedFile> cfile(new CCompressedFile());
	unique_ptr<CVCFIO> vcf_io(new CVCFIO());
	bool end_of_processing = false;

	if (!vcf->OpenForReading(params.vcf_file_name))
	{
		cerr << "Cannot open: " << params.vcf_file_name << endl;
		return false;
	}

	vcf_io->Connect(vcf.get());

	string header;
	vector<string> v_samples;
    
    int no_flt_keys, no_info_keys, no_fmt_keys;

    std::vector<int> FilterIdToFieldId;
    std::vector<int> InfoIdToFieldId;
    std::vector<int> FormatIdToFieldId;

	int gt_key_id = -1; //default: no gt in keys
    vcf->GetFilterInfoFormatKeys(no_flt_keys, no_info_keys, no_fmt_keys,keys, gt_key_id); 
	cfile->SetGTId(gt_key_id);

	cfile->SetNeglectLimit(params.neglect_limit);
	cfile->SetNoSamples(vcf->GetNoSamples());
	cfile->SetPloidy(vcf->GetPloidy());
	cfile->SetNoThreads(params.no_threads);

	size_t i_variant = 0;

    for(int i = 0; i < no_flt_keys+no_info_keys+no_fmt_keys; i++)
    {
        switch(keys[i].keys_type)
        {
        case key_type_t::flt:
			if (keys[i].key_id >= FilterIdToFieldId.size())
				FilterIdToFieldId.resize(keys[i].key_id + 1, -1);
			FilterIdToFieldId[keys[i].key_id] = i;
            break;
        case key_type_t::info:
			if (keys[i].key_id >= InfoIdToFieldId.size())
				InfoIdToFieldId.resize(keys[i].key_id + 1, -1);
            InfoIdToFieldId[keys[i].key_id] = i;
            break;
        case key_type_t::fmt:
			if (keys[i].key_id >= FormatIdToFieldId.size())
				FormatIdToFieldId.resize(keys[i].key_id + 1, -1);
            FormatIdToFieldId[keys[i].key_id] = i;
            break;
        }
    }
    
	vcf->GetHeader(header);
	vcf->GetSamplesList(v_samples);
	cfile->SetHeader(header);
	cfile->AddSamples(v_samples);
    cfile->SetNoKeys((uint32_t)keys.size());
    cfile->SetKeys(keys);
	cfile->SetCompressionLevel(params.vcs_compression_level);
    
	function_data_item_t empty_data_map;

	vector<bool> v_empty(no_flt_keys + no_info_keys + no_fmt_keys, true);

	if (!cfile->OpenForWriting(params.db_file_name, no_flt_keys + no_info_keys + no_fmt_keys))
		return false;

	v_bcf_io.reserve(no_variants_in_buf);
	v_bcf_parse.reserve(no_variants_in_buf);
	s_bcf_io = 0;
	s_bcf_parse = 0;

	for (size_t i = 0; i < no_variants_in_buf; ++i)
	{
		v_bcf_io.emplace_back(vcf_io->InitRecord());
		v_bcf_parse.emplace_back(vcf_io->InitRecord());
	}

	unique_ptr<thread> t_io(new thread([&] {
		while (!end_of_processing)
		{
			s_bcf_io = 0;
			for (size_t i = 0; i < no_variants_in_buf; ++i, ++s_bcf_io)
				if (!vcf_io->LoadRecord(v_bcf_io[i]))
					break;

//			cout << "s_bcf_io " + to_string(s_bcf_io) + "\n";

			barrier.count_down_and_wait();
			barrier.count_down_and_wait();
		}
	}));


	// Thread reading VCF files in parts
	unique_ptr<thread> t_vcf(new thread([&] {
		barrier.count_down_and_wait();
		barrier.count_down_and_wait();

		while (!end_of_processing)
		{
//			vector<pair<int, int>> data_to_remove, size_to_remove;
//			auto t1 = high_resolution_clock::now();

			v_vcf_data_io.clear();
//			for (size_t i = 0; i < no_variants_in_buf; ++i)
			for(size_t i = 0; i < s_bcf_parse; ++i)
			{
//				v_vcf_data_io.push_back(make_pair(variant_desc_t(), vector<field_desc>(keys.size())));
				v_vcf_data_io.emplace_back(variant_desc_t(), vector<field_desc>(keys.size()));
//				if (!vcf->GetVariant(v_vcf_data_io.back().first, v_vcf_data_io.back().second, FilterIdToFieldId, InfoIdToFieldId, FormatIdToFieldId))
				if (!vcf->GetVariantFromRec(v_bcf_parse[i], v_vcf_data_io.back().first, v_vcf_data_io.back().second, FilterIdToFieldId, InfoIdToFieldId, FormatIdToFieldId))
				{
					v_vcf_data_io.pop_back();
					break;
				}

#if 0
				// Graph update
				auto& vec = v_vcf_data_io.back().second;

				// Look for empty
				for (int j = 0; j < no_flt_keys + no_info_keys + no_fmt_keys; ++j)
					if (vec[j].present)
						v_empty[j] = false;

				data_to_remove.clear();

				// Not used at this version - for future extension
				for (auto& edge : function_data_graph)
				{
					int from = edge.first.first;
					int to = edge.first.second;

					if (vec[from].data_size > 4 || vec[to].data_size > 4)
					{
						data_to_remove.emplace_back(from, to);
						continue;
					}

					vector<uint8_t> from_val(vec[from].data, vec[from].data + vec[from].data_size * 4);
					vector<uint8_t> to_val(vec[to].data, vec[to].data + vec[to].data_size * 4);

					// !!! Test
					if(from_val != to_val)
					{
						data_to_remove.emplace_back(from, to);
						continue;
					}

					auto p = edge.second.find(from_val);
					if (p == edge.second.end())
						edge.second.emplace(from_val, to_val);
					else
						if (p->second != to_val)
							data_to_remove.emplace_back(from, to);

					if(edge.second.size() > max_size_of_function)
						data_to_remove.emplace_back(from, to);
				}

				for (auto& x : data_to_remove)
					function_data_graph.erase(x);
#endif
			}

/*			auto t2 = high_resolution_clock::now();
			duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

			std::cout << string("    Reader                         ") + to_string(time_span.count()) + "\n";
			*/

			barrier.count_down_and_wait();
			barrier.count_down_and_wait();
        }
	}));

	// Thread making PBWT and compressing data
	unique_ptr<thread> t_compress(new thread([&] {
		barrier.count_down_and_wait();
		barrier.count_down_and_wait();

		while (!end_of_processing)
		{
//			auto t1 = high_resolution_clock::now();
	
			for (size_t i = 0; i < v_vcf_data_compress.size(); ++i)
            {
				i_variant++;

				cfile->SetVariant(v_vcf_data_compress[i].first, v_vcf_data_compress[i].second);

                for(size_t j = 0; j < keys.size(); ++j)
                {
                    if(v_vcf_data_compress[i].second[j].data_size > 0)
                    {
						if(v_vcf_data_compress[i].second[j].data)
							delete[] v_vcf_data_compress[i].second[j].data;
                        v_vcf_data_compress[i].second[j].data = nullptr;
                        v_vcf_data_compress[i].second[j].data_size = 0;
                    }
                }
            }

/*			auto t2 = high_resolution_clock::now();
			duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

			std::cout << string("*** Storer  ") + to_string(time_span.count()) + "\n";*/

            barrier.count_down_and_wait();
			barrier.count_down_and_wait();
		}
	}));

	// Synchronization
	size_t no_variants = 0;
	while (!end_of_processing)
	{
		barrier.count_down_and_wait();
		
		swap(v_bcf_io, v_bcf_parse);
		s_bcf_parse = s_bcf_io;

		swap(v_vcf_data_compress, v_vcf_data_io);
		if (v_vcf_data_compress.empty() && s_bcf_parse == 0)
			end_of_processing = true;

		barrier.count_down_and_wait();

		no_variants += v_vcf_data_compress.size();
		cout << no_variants << "\r";
//		cout << "No variants" + to_string(no_variants) + "\n";
		fflush(stdout);
	}

	t_compress->join();
	t_vcf->join();
	t_io->join();

	for (auto p : v_bcf_io)
		vcf_io->ReleaseRecord(p);
	for (auto p : v_bcf_parse)
		vcf_io->ReleaseRecord(p);

	cfile->Close();

	vcf->Close();
	cout << endl;

	cfile->OptimizeDB(function_size_graph, function_data_graph);

	cout << endl;

	return true;
}

// ******************************************************************************
bool CApplication::DecompressDB()
{
	CBarrier barrier(4);
	unique_ptr<CVCF> vcf(new CVCF());
	unique_ptr<CCompressedFile> cfile(new CCompressedFile());
	unique_ptr<CVCFIO> vcf_io(new CVCFIO());
	bool end_of_processing = false;

	if (!vcf->OpenForWriting(params.vcf_file_name, params.out_type, params.bcf_compression_level))
	{
		cerr << "Cannot open: " << params.vcf_file_name << endl;
		return false;
	}

	cfile->SetNoThreads(params.no_threads);

	if (!cfile->OpenForReading(params.db_file_name))
		return false;

	params.neglect_limit = cfile->GetNeglectLimit();

	uint32_t no_variants = cfile->GetNoVariants();
	uint32_t i_variant = 0;

	string header;
	vector<string> v_samples;

	cfile->GetHeader(header);
	cfile->GetSamples(v_samples);
    cfile->GetKeys(keys);
	vcf->SetHeader(header);
	vcf->AddSamples(v_samples);
	vcf->WriteHeader();
	vcf->SetPloidy(cfile->GetPloidy());

	vcf_io->Connect(vcf.get());

	v_bcf_io.reserve(no_variants_in_buf);
	v_bcf_parse.reserve(no_variants_in_buf);
	s_bcf_io = 0;
	s_bcf_parse = 0;

	for (size_t i = 0; i < no_variants_in_buf; ++i)
	{
		v_bcf_io.emplace_back(vcf_io->InitRecord());
		v_bcf_parse.emplace_back(vcf_io->InitRecord());
	}

	// Thread for low level I/O for VCF file
	unique_ptr<thread> t_io(new thread([&] {
		barrier.count_down_and_wait();
		barrier.count_down_and_wait();

		while (!end_of_processing)
		{
			for (size_t i = 0; i < s_bcf_io; ++i)
				if (!vcf_io->StoreRecord(v_bcf_io[i]))
					break;

			//			cout << "s_bcf_io " + to_string(s_bcf_io) + "\n";

			barrier.count_down_and_wait();
			barrier.count_down_and_wait();
		}
		}));

	// Thread making rev-PBWT and decompressing data
	unique_ptr<thread> t_compress(new thread([&] {
		while (!end_of_processing)
		{
			v_vcf_data_compress.clear();

			for (size_t i = 0; i < no_variants_in_buf && i_variant < no_variants; ++i, ++i_variant)
			{
//				v_vcf_data_compress.push_back(make_pair(variant_desc_t(), vector<field_desc>(keys.size())));
				v_vcf_data_compress.emplace_back(variant_desc_t(), vector<field_desc>(keys.size()));
				cfile->GetVariant(v_vcf_data_compress.back().first, v_vcf_data_compress.back().second);
			}
			
			barrier.count_down_and_wait();
			barrier.count_down_and_wait();
		}
	}));

	// Thread writing VCF files in parts
	unique_ptr<thread> t_vcf(new thread([&] {
		while (!end_of_processing)
		{
			s_bcf_parse = 0;
			for (size_t i = 0; i < v_vcf_data_io.size(); ++i, ++s_bcf_parse)
            {
                vcf->SetVariantToRec(v_bcf_parse[i], v_vcf_data_io[i].first, v_vcf_data_io[i].second, keys);
                
                for(size_t j = 0; j < keys.size(); ++j)
                    if(v_vcf_data_io[i].second[j].data_size)
                    {
                        delete[] v_vcf_data_io[i].second[j].data;
                        v_vcf_data_io[i].second[j].data = nullptr;
                        v_vcf_data_io[i].second[j].data_size = 0;
                    }
            }
            v_vcf_data_io.clear();

			barrier.count_down_and_wait();
			barrier.count_down_and_wait();
		}
	}));

	// Synchronization
	while (!end_of_processing)
	{
		barrier.count_down_and_wait();

		swap(v_bcf_io, v_bcf_parse);
		s_bcf_io = s_bcf_parse;

		swap(v_vcf_data_compress, v_vcf_data_io);
		if (v_vcf_data_io.empty() && s_bcf_io == 0)
			end_of_processing = true;

		cout << i_variant << "\r";
		fflush(stdout);
		barrier.count_down_and_wait();
	}

	t_compress->join();
	t_vcf->join();
	t_io->join();

	for (auto p : v_bcf_io)
		vcf_io->ReleaseRecord(p);
	for (auto p : v_bcf_parse)
		vcf_io->ReleaseRecord(p);

	cfile->Close();
	vcf->Close();
	cout << endl;

	return true;
}

// EOF
