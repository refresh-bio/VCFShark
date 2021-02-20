#pragma once
// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <list>
#include <deque>
#include <memory>

#include "params.h"
#include "vcf.h"
#include "cfile.h"
#include "vcf.h"

using namespace std;

// ******************************************************************************
class CApplication
{
	const size_t no_variants_in_buf = 8192u;
//	const size_t no_variants_in_buf = 4096u;
	const size_t max_size_of_function = 16384u;

	typedef pair<uint8_t, uint32_t> run_desc_t;

	list<vector<run_desc_t>> l_hist_rle_genotypes;
	list<array<uint8_t, 2>> l_hist_tracked_samples;

	deque<vector<run_desc_t>> d_hist_rle_genotypes;
	deque<array<uint8_t, 2>> d_hist_tracked_samples;

	CParams params;

	vector<pair<variant_desc_t, vector<field_desc>>> v_vcf_data_compress, v_vcf_data_io;

	vector<tuple<uint8_t, run_t, uint32_t, uint32_t>> v_sample_data_compress, v_sample_data_io;
	vector<pair<variant_desc_t, uint8_t>> v_sample_d_data_compress, v_sample_d_data_io;

	vector<bcf1_t*> v_bcf_io, v_bcf_parse;
	size_t s_bcf_io, s_bcf_parse;

	function_size_graph_t function_size_graph;
	function_data_graph_t function_data_graph;

	mutex mtx;
	condition_variable cv;

    vector<key_desc> keys;

public:
	CApplication(const CParams &_params);
	~CApplication();

	bool CompressDB();
	bool DecompressDB();
};

// EOF
