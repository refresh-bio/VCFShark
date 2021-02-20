#pragma once
// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include "defs.h"
#include <mutex>
#include <condition_variable>
#include <nmmintrin.h>
#include <string>
#include <vector>

using namespace std;

#ifdef OUR_STRTOL
#ifdef __APPLE__
long int strtol(const char* str, char** endptr, int base) ;
#else
long int strtol(const char* str, char** endptr, int base) noexcept;
#endif
#endif

// *****************************************************************************************
void cumulate_sums(array<uint32_t, SIGMA> &a_hist, uint32_t &max_count);
void cumulate_sums(vector<uint32_t> &v_hist, uint32_t &max_count);


// ************************************************************************************
template <typename T>
void calc_cumulate_histogram(const vector<T>& data, vector<uint32_t>& v_hist, uint32_t& max_count)
{
	fill(v_hist.begin(), v_hist.end(), 0u);

	for (auto x : data)
		++v_hist[x];

	cumulate_sums(v_hist, max_count);
}

// ************************************************************************************
template <typename T>
void calc_cumulate_histogram(const vector<pair<T, uint32_t>>& rle_data, vector<uint32_t>& v_hist, uint32_t& max_count)
{
	fill(v_hist.begin(), v_hist.end(), 0u);

	for (auto x : rle_data)
		v_hist[x.first] += x.second;

	cumulate_sums(v_hist, max_count);
}

// *****************************************************************************************
//
class CBarrier
{
public:
	CBarrier(const CBarrier&) = delete;
	CBarrier& operator=(const CBarrier&) = delete;
	explicit CBarrier(unsigned int count) :
		m_count(count), m_generation(0),
		m_count_reset_value(count)
	{
	}
	void count_down_and_wait()
	{
		std::unique_lock< std::mutex > lock(m_mutex);
		unsigned int gen = m_generation;
		if (--m_count == 0)
		{
			m_generation++;
			m_count = m_count_reset_value;
			m_cond.notify_all();
			return;
		}
		while (gen == m_generation)
			m_cond.wait(lock);
	}
private:
	std::mutex m_mutex;
	std::condition_variable m_cond;
	unsigned int m_count;
	unsigned int m_generation;
	unsigned int m_count_reset_value;
};

// *****************************************************************************************
template<typename T> T NormalizeValue(T val, T min_val, T max_val)
{
	if (val < min_val)
		return min_val;
	if (val > max_val)
		return max_val;
	return val;
}

// *****************************************************************************************
template<typename T>
uint32_t pop_count(T x)
{
	uint32_t r = 0;

	for (; x; ++r)
		x &= x - 1;

	return r;
}

// *****************************************************************************************
template<typename T>
uint32_t ilog2(T x)
{
	uint32_t r = 0;

	for (; x; ++r)
		x >>= 1;

	return r;
}

// *****************************************************************************************
template<typename T>
uint64_t no_bytes(T x)
{
	uint64_t r = 1;

	x >>= 8;

	for (; x; ++r)
		x >>= 8;

	return r;
}

uint64_t popcnt(uint64_t x);
string trim(string s);

// *****************************************************************************************
template <typename T>
uint64_t modulo_divisor(T x, int mod)
{
	switch (mod)
	{
	case 0: return 0;
	case 1: return 0;
	case 2: return x % 2;
	case 3: return x % 3;
	case 4: return x % 4;
	case 5: return x % 5;
	case 6: return x % 6;
	case 7: return x % 7;
	case 8: return x % 8;
	case 9: return x % 9;
	case 10: return x % 10;
	case 11: return x % 11;
	case 12: return x % 12;
	case 13: return x % 13;
	case 14: return x % 14;
	case 15: return x % 15;
	case 16: return x % 16;
	case 17: return x % 17;
	case 18: return x % 18;
	case 19: return x % 19;
	case 20: return x % 20;
	case 21: return x % 21;
	case 22: return x % 22;
	case 23: return x % 23;
	case 24: return x % 24;
	case 25: return x % 25;
	case 26: return x % 26;
	case 27: return x % 27;
	case 28: return x % 28;
	case 29: return x % 29;
	case 30: return x % 30;
	case 31: return x % 31;
	case 32: return x % 32;
	case 33: return x % 33;
	case 34: return x % 34;
	case 35: return x % 35;
	case 36: return x % 36;
	case 37: return x % 37;
	case 38: return x % 38;
	case 39: return x % 39;
	case 40: return x % 40;
	case 41: return x % 41;
	case 42: return x % 42;
	case 43: return x % 43;
	case 44: return x % 44;
	case 45: return x % 45;
	case 46: return x % 46;
	case 47: return x % 47;
	case 48: return x % 48;
	case 49: return x % 49;
	case 50: return x % 50;
	case 51: return x % 51;
	case 52: return x % 52;
	case 53: return x % 53;
	case 54: return x % 54;
	case 55: return x % 55;
	case 56: return x % 56;
	case 57: return x % 57;
	case 58: return x % 58;
	case 59: return x % 59;
	case 60: return x % 60;
	case 61: return x % 61;
	case 62: return x % 62;
	case 63: return x % 63;
	case 64: return x % 64;
	default: return x % mod;
	}
}

// EOF
