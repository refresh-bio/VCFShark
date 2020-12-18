#pragma once
// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2020-12-18
// *******************************************************************************************

#include "params.h"
#include <vector>
#include <string>
#include <unordered_map>
#include <htslib/hts.h>
#include <htslib/vcf.h>

using namespace std;

enum  class key_type_t {flt, info, fmt};  // FILTER / INFO / FORMAT

// ************************************************************************************
typedef struct key_desc_tag {
    uint32_t key_id;
    key_type_t keys_type;
    int8_t type; //one of: BCF_HT_FLAG 0 / BCF_HT_INT  1 / BCF_HT_REAL 2 / BCF_HT_STR 3
} key_desc;


// ************************************************************************************
typedef struct field_desc_tag {
    bool present = false;  // true if present in description
    char * data = nullptr;
    uint32_t data_size = 0; //current size of allocated memory

} field_desc;

typedef struct variant_desc_tag {
	string chrom;
	int64_t pos;
	string id;
	string ref;
	string alt;
	string qual;
	uint32_t n_allele;

	bool operator==(const struct variant_desc_tag &x)
	{
		if (chrom == "" && x.chrom == "")
			return true;

		return chrom == x.chrom &&
			pos == x.pos;
	}

	bool operator !=(const struct variant_desc_tag &x)
	{
		return !operator==(x);
	}

	bool operator<(const struct variant_desc_tag &x)
	{
		if (chrom != x.chrom)
		{
			if (chrom.empty())
				return false;
			if (x.chrom.empty())
				return true;

			return chrom < x.chrom;
		}
		if (pos != x.pos)
			return pos < x.pos;
		return alt < x.alt;
	}
} variant_desc_t;

// ************************************************************************************
class CVCF
{
    htsFile * vcf_file;
    bcf1_t * rec;
    int ploidy;
    
    bool first_variant;
    int curr_alt_number; //allele from ALT field (from 1)

    void *dst_int = nullptr;
    void *dst_real = nullptr;
    void *dst_str = nullptr;
    void *dst_flag = nullptr;
    int  ndst_int = 0;
    int  ndst_real = 0;
    int  ndst_str = 0;
    int  ndst_flag = 0;
    
public:
    
    bcf_hdr_t * vcf_hdr;
    
	CVCF();
	~CVCF();

	// Open VCF file for reading
	bool OpenForReading(string & file_name);

	// Open VCF file for writing
	bool OpenForWriting(string & file_name, file_type type, char bcf_compression_level);

	// Close VCF file
	bool Close();

	// If open, return no. of samples
	int GetNoSamples();

	// Get complete header as a string
	bool GetHeader(string &v_header);
	
	// Parse complete header given as a string
	bool SetHeader(string & v_header);
	
	// If open for writing store the header
	bool WriteHeader();
	
	// Get ploidy
	int GetPloidy();
    
	// Set ploidy
	void SetPloidy(int _ploidy);
    
    // If open, return no. of possible FLT/INFO/FORMAT fields and the keys in vector keys
    bool GetFilterInfoFormatKeys(int &no_flt_keys, int &no_info_keys, int &no_fmt_keys, vector<key_desc> &keys, int & gt_key_id);
    
	// If file open give the next variant:
	// desc - variant description
	// data - genotypes encoded (1B for each sample) as:
	// * bits 0-1 information about 1st haplotype:
	//     00 - absent
	//     01 - present
	//     10 - multi-allele
	//     11 - unknown
	// * bits 2-3 - information about 2nd haplotype, encoded as above
	// * bit 4:
	//     0 - no phasing info
	//     1 - data phase
	bool GetVariant(variant_desc_t &desc, vector<field_desc> &fields,  std::unordered_map<int, uint32_t> &FilterIdToFieldId, std::unordered_map<int, uint32_t> &InfoIdToFieldId, 
		std::unordered_map<int, uint32_t> &FormatIdToFieldId);

	// Store info about variant - parameters the same as for GetVariant
	bool SetVariant(variant_desc_t &desc, vector<field_desc> &fields, vector<key_desc> keys);
	
	// Get vector with sample names
	bool GetSamplesList(vector<string> &s_list);
	
	// Add sample names to file
	bool AddSamples(vector<string> &s_list);
    
	// Add a single sample name
	bool AddSample(string &s_name);
};

// EOF
