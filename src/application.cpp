// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2020-12-18
// *******************************************************************************************

#include "application.h"
#include "utils.h"
#include "graph_opt.h"

#include <iostream>
#include <vector>

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
	CBarrier barrier(3);
	unique_ptr<CVCF> vcf(new CVCF());
	unique_ptr<CCompressedFile> cfile(new CCompressedFile());
	bool end_of_processing = false;

	if (!vcf->OpenForReading(params.vcf_file_name))
	{
		cerr << "Cannot open: " << params.vcf_file_name << endl;
		return false;
	}

	string header;
	vector<string> v_samples;
    
    int no_flt_keys, no_info_keys, no_fmt_keys;

    std::unordered_map<int, uint32_t> FilterIdToFieldId;
    std::unordered_map<int, uint32_t> InfoIdToFieldId;
    std::unordered_map<int, uint32_t> FormatIdToFieldId;
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
			FilterIdToFieldId[keys[i].key_id] = i;
            break;
        case key_type_t::info:
            InfoIdToFieldId[keys[i].key_id] = i;
            break;
        case key_type_t::fmt:
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
    
	function_data_item_t empty_data_map;

	vector<bool> v_empty(no_flt_keys + no_info_keys + no_fmt_keys, true);

	if (!cfile->OpenForWriting(params.db_file_name, no_flt_keys + no_info_keys + no_fmt_keys))
		return false;

	// Thread reading VCF files in parts
	unique_ptr<thread> t_io(new thread([&] {
		while (!end_of_processing)
		{
			vector<pair<int, int>> data_to_remove, size_to_remove;

			v_vcf_data_io.clear();
			for (size_t i = 0; i < no_variants_in_buf; ++i)
			{
				v_vcf_data_io.push_back(make_pair(variant_desc_t(), vector<field_desc>(keys.size())));
				if (!vcf->GetVariant(v_vcf_data_io.back().first, v_vcf_data_io.back().second, FilterIdToFieldId, InfoIdToFieldId, FormatIdToFieldId))
				{
					v_vcf_data_io.pop_back();
					break;
				}

				// Aktualizacja grafu
				auto& vec = v_vcf_data_io.back().second;

				// Wyszukiwanie pustych
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

			}

			barrier.count_down_and_wait();
			barrier.count_down_and_wait();
        }
	}));

	// Thread making PBWT and compressing data
	unique_ptr<thread> t_vcf(new thread([&] {
		while (!end_of_processing)
		{
			for (size_t i = 0; i < v_vcf_data_compress.size(); ++i)
            {
				i_variant++;

				cfile->SetVariant(v_vcf_data_compress[i].first, v_vcf_data_compress[i].second);

                for(size_t j = 0; j < keys.size(); ++j)
                {
                    if(v_vcf_data_compress[i].second[j].data_size > 0)
                    {
                        delete[] v_vcf_data_compress[i].second[j].data;
                        v_vcf_data_compress[i].second[j].data = nullptr;
                        v_vcf_data_compress[i].second[j].data_size = 0;
                    }
                }
            }
            barrier.count_down_and_wait();
			barrier.count_down_and_wait();
		}
	}));

	// Synchronization
	size_t no_variants = 0;
	while (!end_of_processing)
	{
		barrier.count_down_and_wait();
		swap(v_vcf_data_compress, v_vcf_data_io);
		if (v_vcf_data_compress.empty())
			end_of_processing = true;

		no_variants += v_vcf_data_compress.size();
		cout << no_variants << "\r";
		fflush(stdout);
		barrier.count_down_and_wait();
	}

	t_vcf->join();
	t_io->join();

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
	CBarrier barrier(3);
	unique_ptr<CVCF> vcf(new CVCF());
	unique_ptr<CCompressedFile> cfile(new CCompressedFile());
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

	// Thread making rev-PBWT and decompressing data
	unique_ptr<thread> t_vcf(new thread([&] {
		while (!end_of_processing)
		{
			v_vcf_data_compress.clear();

			for (size_t i = 0; i < no_variants_in_buf && i_variant < no_variants; ++i, ++i_variant)
			{
				v_vcf_data_compress.push_back(make_pair(variant_desc_t(), vector<field_desc>(keys.size())));
				cfile->GetVariant(v_vcf_data_compress.back().first, v_vcf_data_compress.back().second);
			}
			
			barrier.count_down_and_wait();
			barrier.count_down_and_wait();
		}
	}));

	// Thread writing VCF files in parts
	unique_ptr<thread> t_io(new thread([&] {
		while (!end_of_processing)
		{
			for (size_t i = 0; i < v_vcf_data_io.size(); ++i)
            {
                vcf->SetVariant(v_vcf_data_io[i].first, v_vcf_data_io[i].second, keys);
                
                for(size_t j = 0; j< keys.size(); ++j)
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
		swap(v_vcf_data_compress, v_vcf_data_io);
		if (v_vcf_data_io.empty())
			end_of_processing = true;

		cout << i_variant << "\r";
		fflush(stdout);
		barrier.count_down_and_wait();
	}

	t_vcf->join();
	t_io->join();

	cfile->Close();
	vcf->Close();
	cout << endl;

	return true;
}

// EOF
