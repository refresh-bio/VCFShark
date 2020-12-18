#pragma once
// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2020-12-18
// *******************************************************************************************

#include <vector>
#include <cstdint>
#include "defs.h"

using namespace std;

// ************************************************************************************
class CBuffer {
public:
	enum class buffer_t {none, flag, integer, real, text};

private:
	uint32_t max_size;
	buffer_t type;

	vector<uint32_t> v_size;
	vector<uint8_t> v_data;

	bool is_function;
	bool is_no_data;		// buffer was set to empty vector

	function_data_item_t fun;

	uint32_t v_size_pos;
	uint32_t v_data_pos;

	uint32_t encode_var_int(char* p);
	uint32_t decode_var_int(uint8_t* p, uint32_t &val);
	uint32_t encode_float(char* p);
	uint32_t decode_float(uint8_t* p, float &val);
	void permute_integer_series_forward();
	void permute_integer_series_backward();
	void permute_float_series_forward();
	void permute_float_series_backward();

	int gcd(int n, int m);

public:
	CBuffer();
	~CBuffer();

	void SetMaxSize(uint32_t _max_size);

	// Output buffer methods
	void WriteFlag(uint8_t f);
	void WriteInt(char* p, uint32_t size);
	void WriteIntVarSize(char* p, uint32_t size);
	void WriteInt64(int64_t x);
	void WriteReal(char* p, uint32_t size);
	void WriteText(char* p, uint32_t size);
	void GetBuffer(vector<uint32_t>& _v_size, vector<uint8_t>& _v_data);
	bool IsFull(void);

	// Input buffer methods
	void ReadFlag(uint8_t &flag);
	void ReadInt(char* &p, uint32_t& size);
	void ReadIntVarSize(char* &p, uint32_t& size);
	void ReadInt64(int64_t &x);
	void ReadReal(char* &p, uint32_t& size);
	void ReadText(char* &p, uint32_t& size);
	void SetBuffer(vector<uint32_t>& _v_size, vector<uint8_t>& _v_data);

	void SetFunction(function_data_item_t& _fun);
	void FuncInt(char*& p, uint32_t& size, char* src_p, uint32_t src_size);
	void FuncReal(char*& p, uint32_t& size, char* src_p, uint32_t src_size);

	bool IsEmpty();
};

// EOF
