// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2020-12-18
// *******************************************************************************************

#include "vcf.h"
#include <iostream>
#include <cassert>

// ************************************************************************************
CVCF::CVCF()
{
    vcf_file = nullptr;
    vcf_hdr = nullptr;
    rec = nullptr;
    curr_alt_number = 1;
    ploidy = 0; //default
    first_variant = true;
}

// ************************************************************************************
CVCF::~CVCF()
{
    Close();
}

// ************************************************************************************
bool CVCF::OpenForReading(string & file_name)
{
    if(vcf_file)
        hts_close(vcf_file);
    vcf_file = hts_open(file_name.c_str(), "r"); //  With 'r' opens for reading; any further format mode letters are ignored as the format is detected by checking the first few bytes or BGZF blocks of the file.
    if(!vcf_file)
        return false;
    hts_set_opt(vcf_file, HTS_OPT_CACHE_SIZE, 32000000);
    if(vcf_hdr)
        bcf_hdr_destroy(vcf_hdr);
    vcf_hdr = bcf_hdr_read(vcf_file);
    rec = bcf_init();
    return true;
}

// ************************************************************************************
bool CVCF::OpenForWriting(string & file_name, file_type type, char bcf_compression_level)
{
    if(type == file_type::VCF)
        vcf_file = hts_open(file_name.c_str(), "w");
    else  // file_type::BCF
    {
        char write_mode[5] = "wb";
        write_mode[2] = bcf_compression_level;
        write_mode[3] = '\0';
        vcf_file = hts_open(file_name.c_str(), write_mode);
    }
    if(!vcf_file)
        return false;
    hts_set_opt(vcf_file, HTS_OPT_CACHE_SIZE, 32000000);
    rec = bcf_init();
    return true;
}

// ************************************************************************************
bool CVCF::Close()
{
    if (vcf_file)
    {
        if (hts_close(vcf_file) < 0)
        {
            cerr << "Cannot close VCF file\n";
            return false;
        }
        vcf_file = nullptr;
    }

    if(vcf_hdr)
    {
        bcf_hdr_destroy(vcf_hdr);
        vcf_hdr = nullptr;
    }

    if(rec)
    {
        bcf_clear(rec);
        bcf_destroy(rec);
        rec = nullptr;
    }

    if (dst_int)
    {
        free(dst_int);
        dst_int = nullptr;
    }

    if (dst_real)
    {
        free(dst_real);
        dst_real = nullptr;
    }

    if (dst_str)
    {
        free(dst_str);
        dst_str = nullptr;
    }

    if (dst_flag)
    {
        free(dst_flag);
        dst_flag = nullptr;
    }

    return true;
}

// ************************************************************************************
int CVCF::GetNoSamples()
{
    if(!vcf_file || !vcf_hdr)
        return -1;
    
    return bcf_hdr_nsamples(vcf_hdr);
}

// ************************************************************************************
bool CVCF::GetSamplesList(vector<string> &s_list)
{
    if(!vcf_hdr)
        return false;
    int n = GetNoSamples();
    for (int i = 0; i < n; i++)
        s_list.push_back(vcf_hdr->samples[i]);

    return true;
}

// ************************************************************************************
int CVCF::GetPloidy()
{
    return ploidy;
}

// ************************************************************************************
void CVCF::SetPloidy(int _ploidy)
{
    ploidy = _ploidy;
}

// ************************************************************************************
bool CVCF::GetHeader(string &v_header)
{
    if(!vcf_hdr)
        return false;
    kstring_t str = {0, 0, nullptr};
    bcf_hdr_format(vcf_hdr, 0, &str);
    char * ptr = strstr(str.s, "#CHROM");
    v_header.assign(str.s, ptr - str.s);
    if(str.m)
        free(str.s);
    return true;
}

// ************************************************************************************
bool CVCF::SetHeader(string &v_header)
{
    vcf_hdr = bcf_hdr_init("r");
    string temp = v_header + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    bcf_hdr_parse(vcf_hdr, (char *)temp.c_str());
    if(bcf_hdr_sync(vcf_hdr) != 0)
        return false;
    if(!vcf_hdr)
        return false;
    return true;
}

// ************************************************************************************
bool CVCF::AddSamples(vector<string> &s_list)
{
    if(!vcf_hdr)
        return false;
    for(size_t i = 0; i < s_list.size(); i++)
        bcf_hdr_add_sample(vcf_hdr, s_list[i].c_str());
    
    if(bcf_hdr_sync(vcf_hdr) != 0)
        return false;
    return true;
}

// ************************************************************************************
bool CVCF::AddSample(string & s_name)
{
    if(!vcf_hdr)
        return false;
    bcf_hdr_add_sample(vcf_hdr, s_name.c_str());
    if(bcf_hdr_sync(vcf_hdr) != 0)
        return false;
    return true;
}

// ************************************************************************************
bool CVCF::WriteHeader()
{
    if(vcf_hdr && vcf_file)
    {
        if(bcf_hdr_write(vcf_file, vcf_hdr) < 0)
            return false;
        return true;
    }
    return false;
}

bool CVCF::GetFilterInfoFormatKeys(int &no_flt_keys, int &no_info_keys, int &no_fmt_keys, vector<key_desc> &keys, int & gt_key_id)
{
    if(!vcf_file || !vcf_hdr)
        return false;
    
    no_flt_keys = 0;
    no_info_keys = 0;
    no_fmt_keys = 0;
    
    key_desc new_key;
    for(int i = 0; i < vcf_hdr->nhrec; i++)
    {
        
        
        if(vcf_hdr->hrec[i]->type == BCF_HL_FLT || vcf_hdr->hrec[i]->type ==  BCF_HL_INFO || vcf_hdr->hrec[i]->type == BCF_HL_FMT)
        {
            //checking if it is a duplicate; if so, curr_dict_id different (and not increased)
            int id = bcf_hdr_id2int(vcf_hdr, BCF_DT_ID, vcf_hdr->hrec[i]->vals[0]);
            new_key.key_id = id;

            if(vcf_hdr->hrec[i]->type == BCF_HL_FLT)
            {
                no_flt_keys++;
                new_key.keys_type = key_type_t::flt;
                new_key.type = BCF_HT_FLAG;
            }
            else if(vcf_hdr->hrec[i]->type ==  BCF_HL_INFO || vcf_hdr->hrec[i]->type == BCF_HL_FMT)
            {
                if(vcf_hdr->hrec[i]->type == BCF_HL_FMT) //FORMAT
                {
                    no_fmt_keys++;
                    new_key.keys_type = key_type_t::fmt;
                    new_key.type = bcf_hdr_id2type(vcf_hdr,BCF_HL_FMT,id);
                     if (strcmp(vcf_hdr->id[BCF_DT_ID][id].key, "GT") == 0)
                     {
                         new_key.type = BCF_HT_INT; // GT is treated as an int
                         gt_key_id = (int) keys.size(); 
                    }
                }
                else //INFO
                {
                    no_info_keys++;
                    new_key.keys_type = key_type_t::info;
                    new_key.type = bcf_hdr_id2type(vcf_hdr,BCF_HL_INFO, id);
                }
            }
            keys.push_back(new_key);
        }
    }
    return true;
}
// ************************************************************************************
bool CVCF::GetVariant(variant_desc_t &desc,  vector<field_desc> &fields, std::unordered_map<int, uint32_t> &FilterIdToFieldId, std::unordered_map<int, uint32_t> &InfoIdToFieldId,  
    std::unordered_map<int, uint32_t> &FormatIdToFieldId)
{
	desc.chrom.clear();

	if (!vcf_file || !vcf_hdr)
		return false;

    if(curr_alt_number == 1) // read new line/record
    {
        bcf_clear(rec);
        if(bcf_read(vcf_file, vcf_hdr, rec) == -1)
            return false;
        
        if(rec->errcode)
        {
            std::cerr << "Error in VCF file\n";
            exit(1);
        }
        bcf_unpack((bcf1_t*)rec, BCF_UN_ALL);
    }
    
    desc.chrom = vcf_hdr->id[BCF_DT_CTG][rec->rid].key; // CHROM
    desc.pos = rec->pos + 1;  // POS
    desc.id = rec->d.id ? rec->d.id : "."; // ID
    desc.n_allele = rec->n_allele;
    
    desc.ref.erase(); //REF
    if (rec->n_allele > 0)
        desc.ref = rec->d.allele[0];
    else
        desc.ref = ".";
    
    desc.alt.erase();  // ALT
    if (rec->n_allele > 1) {
        desc.alt += rec->d.allele[1];
        for(int a = 2; a < rec->n_allele; a++)
        {
            desc.alt += ",";
            desc.alt += rec->d.allele[a];
        }
    }
    else
        desc.alt = ".";
        
    desc.qual.erase(); // QUAL
    if ( bcf_float_is_missing(rec->qual) )
        desc.qual = ".";
    else
    {
        kstring_t s = {0,0,0};
        kputd(rec->qual, &s);
        desc.qual += s.s;
        free(s.s);
    }    
    
    //FILTER
    if (rec->d.n_flt) {
        for (int i = 0; i < rec->d.n_flt; ++i) {
            fields[FilterIdToFieldId[rec->d.flt[i]]].present = true;  //DEV, TMP
        }
    }

    int curr_size = 0;
    
    // INFO
    if (rec->n_info) {
        for (int i = 0; i < rec->n_info; ++i) {
            bcf_info_t *z = &rec->d.info[i];
            fields[InfoIdToFieldId[z->key]].present = true;
            if ( !z->vptr )
                continue;
            if (z->key >= vcf_hdr->n[BCF_DT_ID]) {
                hts_log_error("Invalid BCF, the INFO index is too large");
                return -1;
            }
            if (z->len <= 0) continue;
            
            int type = bcf_hdr_id2type(vcf_hdr,BCF_HL_INFO,z->key);
            int curr_size = 0;
            
            switch (type)
            {
            case BCF_HT_INT:
                curr_size = bcf_get_info_values(vcf_hdr, rec, vcf_hdr->id[BCF_DT_ID][z->key].key, &dst_int, &ndst_int, type);
                break;
            case BCF_HT_REAL:
                curr_size = bcf_get_info_values(vcf_hdr, rec, vcf_hdr->id[BCF_DT_ID][z->key].key, &dst_real, &ndst_real, type);
                break;
            case BCF_HT_STR:
                curr_size = bcf_get_info_values(vcf_hdr, rec, vcf_hdr->id[BCF_DT_ID][z->key].key, &dst_str, &ndst_str, type);
                break;
            case BCF_HT_FLAG:
                curr_size = bcf_get_info_values(vcf_hdr, rec, vcf_hdr->id[BCF_DT_ID][z->key].key, &dst_flag, &ndst_flag, type);
                break;
            }
            
            if(curr_size)
             {
                 switch(type)
                 {
                     case BCF_HT_INT:
                         fields[InfoIdToFieldId[z->key]].data = new char[curr_size * 4];
                         memcpy(fields[InfoIdToFieldId[z->key]].data, (char*)dst_int, curr_size * 4);
                         fields[InfoIdToFieldId[z->key]].data_size = (uint32_t) curr_size;
                         break;
                     case BCF_HT_REAL:
                         fields[InfoIdToFieldId[z->key]].data = new char[curr_size*4];
                         memcpy(fields[InfoIdToFieldId[z->key]].data, (char*)dst_real, curr_size*4);
                         fields[InfoIdToFieldId[z->key]].data_size = curr_size;
                         break;
                     case BCF_HT_STR:
                         fields[InfoIdToFieldId[z->key]].data = new char[curr_size];
                         memcpy(fields[InfoIdToFieldId[z->key]].data, (char*)dst_str, curr_size);
                         fields[InfoIdToFieldId[z->key]].data_size = curr_size;
                         break;
                     case BCF_HT_FLAG:
                         fields[InfoIdToFieldId[z->key]].data = nullptr;
                         fields[InfoIdToFieldId[z->key]].data_size = 0;
                         break;
                 }
             }
            else
            {
                cout << "Error getting variant" << endl;
                return false;
            }
        }
    }

    // FORMAT and individual information
    if (rec->n_sample)
    {
        int i; // , j;
        if ( rec->n_fmt)
        {
            bcf_fmt_t *fmt = rec->d.fmt;
            for (i = 0; i < (int)rec->n_fmt; ++i) {
                if ( !fmt[i].p ) continue;
                if ( fmt[i].id<0 ) //!bcf_hdr_idinfo_exists(h,BCF_HL_FMT,fmt[i].id) )
                {
                    hts_log_error("Invalid BCF, the FORMAT tag id=%d not present in the header", fmt[i].id);
                    abort();
                }
                fields[FormatIdToFieldId[fmt[i].id]].present = true;
    
                int bcf_ht_type = bcf_hdr_id2type(vcf_hdr,BCF_HL_FMT,fmt[i].id); //BCF_HT_INT or BCF_HT_REAL or BCF_HT_FLAG or BCF_HT_STR
                if (strcmp(vcf_hdr->id[BCF_DT_ID][fmt[i].id].key, "GT") == 0)
                {
                    curr_size = bcf_get_format_values(vcf_hdr, rec ,"GT", &dst_int, &ndst_int, BCF_HT_INT);
                }
                else
                {
                    switch (bcf_ht_type)
                    {
                    case BCF_HT_INT:
                        curr_size = bcf_get_format_values(vcf_hdr, rec, vcf_hdr->id[BCF_DT_ID][fmt[i].id].key, &dst_int, &ndst_int, bcf_ht_type);
                        break;
                    case BCF_HT_REAL:
                        curr_size = bcf_get_format_values(vcf_hdr, rec, vcf_hdr->id[BCF_DT_ID][fmt[i].id].key, &dst_real, &ndst_real, bcf_ht_type);
                        break;
                    case BCF_HT_STR:
                        curr_size = bcf_get_format_values(vcf_hdr, rec, vcf_hdr->id[BCF_DT_ID][fmt[i].id].key, &dst_str, &ndst_str, bcf_ht_type);
                        break;
                    case BCF_HT_FLAG:
                        curr_size = bcf_get_format_values(vcf_hdr, rec, vcf_hdr->id[BCF_DT_ID][fmt[i].id].key, &dst_flag, &ndst_flag, bcf_ht_type);
                        break;
                    }
                }
                
                if(curr_size)
                {
                    fields[FormatIdToFieldId[fmt[i].id]].data_size = (uint32_t) curr_size;
                    
                    if (bcf_ht_type == BCF_HT_INT || bcf_ht_type == BCF_HT_REAL || strcmp(vcf_hdr->id[BCF_DT_ID][fmt[i].id].key, "GT") == 0)
                    {
                        fields[FormatIdToFieldId[fmt[i].id]].data = new char[curr_size * 4];

                        if(bcf_ht_type == BCF_HT_INT)
                            memcpy(fields[FormatIdToFieldId[fmt[i].id]].data, dst_int, curr_size * 4);
                        else if (bcf_ht_type == BCF_HT_REAL)
                            memcpy(fields[FormatIdToFieldId[fmt[i].id]].data, dst_real, curr_size * 4);
                        else if (bcf_ht_type == BCF_HT_STR)
                            memcpy(fields[FormatIdToFieldId[fmt[i].id]].data, dst_int, curr_size * 4);  // GTs are ints!
                    }
                    else if (bcf_ht_type == BCF_HT_STR)
                    {
                        fields[FormatIdToFieldId[fmt[i].id]].data = new char[curr_size];
                        memcpy(fields[FormatIdToFieldId[fmt[i].id]].data, dst_str, curr_size);
                    }
                    else
                        assert(0);                    
                }            
            }
        }
    }

    if(first_variant)
        first_variant = false;

    return true;
}

// ************************************************************************************
bool CVCF::SetVariant(variant_desc_t &desc, vector<field_desc> &fields, vector<key_desc> keys)
{
    bcf_clear(rec);
    
    string record;
    record = desc.chrom + "\t0\t" + desc.id + "\t" + desc.ref + "\t" + desc.alt + "\t" + desc.qual +"\t" + "." + "\t" + ".";
    kstring_t s;
    s.s = (char*)record.c_str();
    s.m = record.length();
    s.l = 0;
    vcf_parse(&s, vcf_hdr, rec);
    rec->pos = (int32_t) (desc.pos - 1);
  
    for(size_t j = 0; j < fields.size(); j++){
        if(keys[j].keys_type == key_type_t::flt)
        {
            if(fields[j].present)
                bcf_add_filter(vcf_hdr, rec, keys[j].key_id);
        }
    }
    
    for(size_t j = 0; j < keys.size(); j++){
        if(keys[j].keys_type == key_type_t::info)
        {
            if(fields[j].present){
                bcf_update_info(vcf_hdr, rec, bcf_hdr_int2id(vcf_hdr,BCF_DT_ID,keys[j].key_id), fields[j].data, fields[j].data_size, bcf_hdr_id2type(vcf_hdr,BCF_HL_INFO, keys[j].key_id));
            }
        }
    }
    
    //FORMAT
    for(size_t j = 0; j < keys.size(); j++){
        if(keys[j].keys_type == key_type_t::fmt)
        {
            if(fields[j].present){
                if(strcmp(bcf_hdr_int2id(vcf_hdr,BCF_DT_ID,keys[j].key_id),"GT") == 0)
                        bcf_update_format(vcf_hdr, rec, "GT",fields[j].data, fields[j].data_size, BCF_HT_INT);
                    else
                        bcf_update_format(vcf_hdr, rec, bcf_hdr_int2id(vcf_hdr , BCF_DT_ID, keys[j].key_id), fields[j].data, fields[j].data_size, bcf_hdr_id2type(vcf_hdr,BCF_HL_FMT,keys[j].key_id));
            }
        }
    }
    if(bcf_write(vcf_file, vcf_hdr, rec) != 0)
        return false;
    return true;
}

// EOF
