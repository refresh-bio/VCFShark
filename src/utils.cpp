// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include "utils.h"
#include <iostream>
#include <memory>
#include <sstream>
#include <algorithm>

#ifdef OUR_STRTOL
// *****************************************************************************************
#ifdef __APPLE__
long int strtol(const char* str, char** endptr, int base)
#else
long int strtol(const char* str, char** endptr, int base) noexcept
#endif
{
	if (base != 10)
	{
		std::cerr << "unsuported base " << base << std::endl;
		fflush(stdout);
		exit(1);
	}

	long int val = 0;
	char *p = (char*)str;
	bool is_negative = false;

	if (*p == '-')
	{
		is_negative = true;
		++p;
	}

	while (*p >= '0' && *p <= '9') 
	{
		val = val * 10 + (*p++ - '0');
	}

	if (endptr)
		*endptr = p;

	return is_negative ? -val : val;
}
#endif

// ************************************************************************************
void cumulate_sums(vector<uint32_t> &v_hist, uint32_t &max_count)
{
	// Determine cumulated histogram
	uint32_t v = 0;
	uint32_t cum = 0;
	uint32_t hist_size = (uint32_t) v_hist.size();
	max_count = 0;

	for (uint32_t i = 0; i < hist_size; ++i)
	{
		v = v_hist[i];
		v_hist[i] = cum;
		cum += v;

		if (v > max_count)
			max_count = v;
	}
}

// ************************************************************************************
void cumulate_sums(array<uint32_t, SIGMA> &a_hist, uint32_t &max_count)
{
	// Determine cumulated histogram
	uint32_t v = 0;
	uint32_t cum = 0;
	uint32_t hist_size = (uint32_t)a_hist.size();
	max_count = 0;

	for (uint32_t i = 0; i < hist_size; ++i)
	{
		v = a_hist[i];
		a_hist[i] = cum;
		cum += v;

		if (v > max_count)
			max_count = v;
	}
}

// ************************************************************************************
string trim(string s)
{
	while (!s.empty() & (s.back() == '\n' || s.back() == '\r'))
		s.pop_back();

	return s;
}

// ************************************************************************************
uint64_t popcnt(uint64_t x)
{
	return _mm_popcnt_u64(x);

	return 0;
}

// EOF
