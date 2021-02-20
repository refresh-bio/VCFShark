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
#include "sub_rc.h"
#include <cmath>
#include <algorithm>
#include <numeric>

// *******************************************************************************************
//
// *******************************************************************************************
template<unsigned NO_SYMBOLS, unsigned MAX_LOG_COUNTER, unsigned ADDER>
class CSimpleModel
{
	const uint32_t MAX_TOTAL = 1u << MAX_LOG_COUNTER;
	uint32_t stats[NO_SYMBOLS];
	uint32_t total;

	void rescale()
	{
		while (total >= MAX_TOTAL)
		{
			total = 0;
			for (uint32_t i = 0; i < NO_SYMBOLS; ++i)
			{
				stats[i] = (stats[i] + 1) / 2;
				total += stats[i];
			}
		}
	}

public: 
	CSimpleModel() : total(0)
	{
		fill_n(stats, NO_SYMBOLS, 1u);
		total = NO_SYMBOLS;
	};

	~CSimpleModel()
	{
	};

	CSimpleModel(const CSimpleModel &c) = delete;
	CSimpleModel& operator=(const CSimpleModel&) = delete;

	void Init(int *_init_stats)
	{
		if (_init_stats)
			for (uint32_t i = 0; i < NO_SYMBOLS; ++i)
				stats[i] = _init_stats[i];
		else
			fill_n(stats, NO_SYMBOLS, 1);

		total = accumulate(stats, stats + NO_SYMBOLS, 0u);
		rescale();
	}

	void Init(const CSimpleModel &c)
	{
		copy_n(c.stats, NO_SYMBOLS, stats);
		total = accumulate(stats, stats + NO_SYMBOLS, 0u);
	}

	void GetFreq(int symbol, int &sym_freq, int &left_freq, int &totf)
	{
		left_freq = 0;

		switch (symbol)
		{
			case 4: left_freq += stats[3];
			case 3: left_freq += stats[2];
			case 2: left_freq += stats[1];
			case 1: left_freq += stats[0];
			case 0: break;
			default:
				for (int i = 0; i < symbol; ++i)
					left_freq += stats[i];
		}

		sym_freq = stats[symbol];
		totf = total;
	}

	void Update(int symbol)
	{
		stats[symbol] += ADDER;
		total += ADDER;

		if (total >= MAX_TOTAL)
			rescale();
	}

	int GetSym(int left_freq)
	{
		int t = 0;

		for (uint32_t i = 0; i < NO_SYMBOLS; ++i)
		{
			t += stats[i];
			if (t > left_freq)
				return i;
		}

		return -1;
	}

	uint32_t GetTotal()
	{
		return total;
	}

	uint32_t *GetStats()
	{
		return stats;
	}

	void SetStats(uint32_t *stats_to_set)
	{
		total = 0;
		for (uint32_t i = 0; i < NO_SYMBOLS; ++i)
			total += stats[i] = stats_to_set[i];
	}
};

// *******************************************************************************************
//
// *******************************************************************************************
template <unsigned N_SYMBOLS> class CSimpleModelFixedSize
{
	uint32_t max_total;
	uint32_t stats[N_SYMBOLS];
	uint32_t total;
	uint32_t adder;
	size_t no_updates;

	void rescale()
	{
		while (total >= max_total)
		{
			total = 0;
			for (uint32_t i = 0; i < N_SYMBOLS; ++i)
			{
				stats[i] = (stats[i] + 1) / 2;
				total += stats[i];
			}
		}
	}

public:
	CSimpleModelFixedSize(uint32_t _adder = 1) : max_total(0), total(0), adder(_adder), no_updates(0)
	{
	};

	~CSimpleModelFixedSize()
	{
	};

	CSimpleModelFixedSize(const CSimpleModelFixedSize &c) = delete;
	CSimpleModelFixedSize& operator=(const CSimpleModelFixedSize&) = delete;

	void Init(const int *_init_stats, uint32_t _max_total, uint32_t _adder)
	{
		max_total = _max_total;
		adder = _adder;

		if (_init_stats)
			for (uint32_t i = 0; i < N_SYMBOLS; ++i)
				stats[i] = _init_stats[i];
		else
			fill_n(stats, N_SYMBOLS, 1);

		total = accumulate(stats, stats + N_SYMBOLS, 0u);
		rescale();

		no_updates = 0;
	}

	void Init(const CSimpleModelFixedSize &c)
	{
		max_total = c.max_total;
		adder = c.adder;

		copy_n(c.stats, N_SYMBOLS, stats);
		total = accumulate(stats, stats + N_SYMBOLS, 0u);

		no_updates = 0;
	}

	void GetFreq(int symbol, int &sym_freq, int &left_freq, int &totf)
	{
		left_freq = 0;

		if(N_SYMBOLS >= 5)
			for (int i = 0; i < symbol; ++i)
				left_freq += stats[i];
		else if(N_SYMBOLS > 2)
			switch (symbol)
			{
			case 4: left_freq += stats[3];
			case 3: left_freq += stats[2];
			case 2: left_freq += stats[1];
			case 1: left_freq += stats[0];
			case 0: break;
			}
		else
			switch (symbol)
			{
			case 1: left_freq += stats[0];
			case 0: break;
			}

		sym_freq = stats[symbol];
		totf = total;
	}

	void Update(int symbol)
	{
		stats[symbol] += adder;
		total += adder;

		if (total >= max_total)
			rescale();

		++no_updates;
	}

	int GetSym(int left_freq)
	{
		int t = 0;

		for (uint32_t i = 0; i < N_SYMBOLS; ++i)
		{
			t += stats[i];
			if (t > left_freq)
				return i;
		}

		return -1;
	}

	uint32_t GetTotal()
	{
		return total;
	}

	void Merge(uint32_t *stats_to_merge)
	{
		for (uint32_t i = 0; i < N_SYMBOLS; ++i)
		{
			stats[i] += stats_to_merge[i];
			total += stats_to_merge[i];
		}
	}

	void CompleteMerge()
	{
		rescale();
	}

	uint32_t *GetStats()
	{
		return stats;
	}

	void SetStats(uint32_t *stats_to_set)
	{
		total = 0;
		for (uint32_t i = 0; i < N_SYMBOLS; ++i)
			total += stats[i] = stats_to_set[i];
	}

	void GetLogStats(size_t &_no_updates, vector<float> &_v_freq)
	{
		_no_updates = no_updates;

		_v_freq.resize(N_SYMBOLS, 0u);
		auto sum = accumulate(stats, stats + N_SYMBOLS, 0u);

		if (sum)
			for (uint32_t i = 0; i < N_SYMBOLS; ++i)
				_v_freq[i] = (float)stats[i] / sum;
	}

	size_t GetNoUpdates()
	{
		return no_updates;
	}
};


// *******************************************************************************************
//
// *******************************************************************************************
// max size : 256
// max stat. value: (1 << 24) - 1
template<unsigned NO_SYMBOLS, unsigned MAX_LOG_COUNTER, unsigned ADDER>
class CAdjustableModel
{
	const uint32_t value_mask = 0xffffffu;
	const uint32_t symbol_mask = 0xff000000u;
	const uint32_t symbol_shift = 24;
	const double compact_limit_frac = 0.33;
	const uint16_t compact_limit = std::max(static_cast<uint16_t>(NO_SYMBOLS * compact_limit_frac), (uint16_t)4);
	const uint32_t MAX_TOTAL = 1u << MAX_LOG_COUNTER;

	uint32_t *stats;
//	uint32_t stats_capacity;
//	uint32_t stats_size;
	uint16_t stats_capacity;
	uint16_t stats_size;
	uint32_t total;

	constexpr uint32_t pack_sv(uint32_t s, uint32_t v)
	{
		return (s << symbol_shift) + v;
	}

	void rescale()
	{
		while (total >= MAX_TOTAL)
		{
			total = 0;

			if (stats_capacity)
			{
				for (uint32_t i = 0; i < stats_size; ++i)
				{
					auto &x = stats[i];
					uint32_t v = ((x & value_mask) + 1) / 2;
					x = (x & symbol_mask) + v;
					total += v;
				}

				total += NO_SYMBOLS - (uint32_t) stats_size;
			}
			else
			{
				for (uint32_t i = 0; i < NO_SYMBOLS; ++i)
				{
					auto& x = stats[i];
					x = (x + 1) / 2;
					total += x;
				}
			}
		}
	}

	void resize_stats()
	{
		auto old_stats = stats;

/*		if (stats_capacity < 3)
		{
			stats_capacity *= 2;
		}
		else
		{
			stats_capacity = (uint32_t)(stats_capacity * 1.4);
			if (stats_capacity > compact_limit)
				stats_capacity = compact_limit + 1;
		}*/
		switch (stats_capacity)
		{
		case 1: stats_capacity = 2;	break;
		case 2: stats_capacity = 4;	break;
		case 4: stats_capacity = 8;	break;
		case 8: stats_capacity = 12;	break;
		case 12: stats_capacity = 16;	break;
		case 16: stats_capacity = 20;	break;
		case 20: stats_capacity = 24;	break;
		case 24: stats_capacity = 32;	break;
		case 32: stats_capacity = 48;	break;
		case 48: stats_capacity = 64;	break;
		case 64: stats_capacity = 80;	break;
		case 80: stats_capacity = 96;	break;
		default:
			stats_capacity = (uint16_t)(stats_capacity * 1.4);
		}

		if (stats_capacity > compact_limit)
			stats_capacity = (uint16_t) (compact_limit + (uint16_t) 1);

		stats = new uint32_t[stats_capacity];
		copy_n(old_stats, stats_size, stats);
		delete[] old_stats;
	}

public:
	CAdjustableModel() : total(0)
	{
		stats_capacity = 1;
		stats_size = 0;
		stats = new uint32_t[stats_capacity];
	};

	~CAdjustableModel()
	{
		delete[] stats;
	};

	CAdjustableModel(const CAdjustableModel& c) = delete;
	CAdjustableModel& operator=(const CAdjustableModel&) = delete;

	// Warning: _init_stats are ignored in this model
	void Init(int* _init_stats)
	{
		total = NO_SYMBOLS;
		stats_size = 0;
		stats_capacity = 1;

		if (stats)
			delete[] stats;
		stats = new uint32_t[stats_capacity];
	}

	void Init(const CAdjustableModel& c)
	{
		total = c.total;

		if (stats)
			delete[] stats;

		stats_size = c.stats_size;
		stats_capacity = c.stats_capacity;

		if (stats_capacity)
		{
			stats = new uint32_t[stats_capacity];
			copy_n(c.stats, stats_size, stats);
		}
		else
		{
			stats = new uint32_t[NO_SYMBOLS];
			copy_n(c.stats, NO_SYMBOLS, stats);
		}
	}

	void GetFreq(int symbol, int& sym_freq, int& left_freq, int& totf)
	{
		left_freq = 0;

		if (stats_capacity)
		{
			uint32_t cnt = 0;

			sym_freq = 1;

			for (uint32_t i = 0; i < stats_size; ++i)
			{
				auto x = stats[i];
				uint32_t s = x >> symbol_shift;
				uint32_t v = x & value_mask;

				if ((int) s < symbol)
				{
					left_freq += v;
					++cnt;
				}
				else if ((int) s == symbol)
				{
					sym_freq = v;
					break;
				}
				else
					break;
			}

			left_freq += symbol - cnt;
		}
		else
		{
			switch (symbol)
			{
			case 4: left_freq += stats[3];
			case 3: left_freq += stats[2];
			case 2: left_freq += stats[1];
			case 1: left_freq += stats[0];
			case 0: break;
			default:
				for (int i = 0; i < symbol; ++i)
					left_freq += stats[i];
			}

			sym_freq = stats[symbol];
		}

		totf = total;
	}

	void Update(int symbol)
	{
		if (stats_capacity)
		{
			uint32_t i;
			bool expanded = false;

			for (i = 0; i < stats_size; ++i)
			{
				uint32_t s = stats[i] >> symbol_shift;

				if ((int) s == symbol)
					break;
				else if ((int) s > symbol)
				{
					if (stats_size == stats_capacity)
						resize_stats();

					copy_backward(stats + i, stats + stats_size, stats + stats_size + 1);
					stats[i] = pack_sv(symbol, 1);
					++stats_size;

					expanded = true;
					break;
				}
			}

			if (i == stats_size)
			{
				if(stats_size == stats_capacity)
					resize_stats();

				stats[stats_size] = pack_sv(symbol, 1);
				++stats_size;
				expanded = true;
			}

			stats[i] += ADDER;
			total += ADDER;

			if (expanded && stats_size >= compact_limit)
			{
				auto old_stats = stats;

				stats = new uint32_t[NO_SYMBOLS];
				fill_n(stats, NO_SYMBOLS, 1);
				for (uint32_t i = 0; i < stats_size; ++i)
				{
					auto x = old_stats[i];
					uint32_t s = x >> symbol_shift;
					uint32_t v = x & value_mask;

					stats[s] = v;
				}

				stats_capacity = 0;

				delete[] old_stats;
			}
		}
		else
		{
			stats[symbol] += ADDER;
			total += ADDER;
		}

		if (total >= MAX_TOTAL)
			rescale();
	}

	int GetSym(uint32_t left_freq)
	{
		uint32_t t = 0;

		if (stats_capacity)
		{
			uint32_t cnt = 0;

			for (uint32_t i = 0; i < stats_size; ++i)
			{
				auto x = stats[i];
				uint32_t s = x >> symbol_shift;
				uint32_t v = x & value_mask;

				t += v;

				if (t + (s - cnt) > left_freq)
				{
					if (t + (s - cnt) - v <= left_freq)
						return s;
					return s - (t + (s - cnt) - v - left_freq);
				}

				++cnt;
			}

			return NO_SYMBOLS - (total - left_freq);
		}
		else
		{
			for (uint32_t i = 0; i < NO_SYMBOLS; ++i)
			{
				t += stats[i];
				if (t > left_freq)
					return i;
			}
		}

		return -1;
	}

	uint32_t GetTotal()
	{
		return total;
	}
};


// *******************************************************************************************
//
// *******************************************************************************************
// max size : 256
// max stat. value: (1 << 24) - 1
template<unsigned NO_SYMBOLS, unsigned MAX_LOG_COUNTER, unsigned ADDER>
class CAdjustableModelEmb
{
	const uint32_t value_mask = 0xffffffu;
	const uint32_t symbol_mask = 0xff000000u;
	const uint32_t symbol_shift = 24;
	const double compact_limit_frac = 0.33;
	const uint16_t compact_limit = std::max(static_cast<uint16_t>(NO_SYMBOLS * compact_limit_frac), (uint16_t)4);
	const uint32_t MAX_TOTAL = 1u << MAX_LOG_COUNTER;

	union {
		uint32_t* stats;
		uint32_t emb_stats[2];
	} u_stats;
	uint16_t stats_capacity;
	uint16_t stats_size;
	uint32_t total;

	constexpr uint32_t pack_sv(uint32_t s, uint32_t v)
	{
		return (s << symbol_shift) + v;
	}

	void rescale()
	{
		while (total >= MAX_TOTAL)
		{
			total = 0;

			if (stats_capacity == 2)
			{
				if (stats_size == 2)
				{
					auto& x = u_stats.emb_stats[1];
					uint32_t v = ((x & value_mask) + 1) / 2;
					x = (x & symbol_mask) + v;
					total += v;
				}
				if (stats_size >= 1)
				{
					auto& x = u_stats.emb_stats[0];
					uint32_t v = ((x & value_mask) + 1) / 2;
					x = (x & symbol_mask) + v;
					total += v;
				}

				total += NO_SYMBOLS - (uint32_t)stats_size;
			}
			else if (stats_capacity)
			{
				for (uint32_t i = 0; i < stats_size; ++i)
				{
					auto& x = u_stats.stats[i];
					uint32_t v = ((x & value_mask) + 1) / 2;
					x = (x & symbol_mask) + v;
					total += v;
				}

				total += NO_SYMBOLS - (uint32_t)stats_size;
			}
			else
			{
				for (uint32_t i = 0; i < NO_SYMBOLS; ++i)
				{
					auto& x = u_stats.stats[i];
					x = (x + 1) / 2;
					total += x;
				}
			}
		}
	}

	void resize_stats()
	{
		if (stats_capacity == 2)
		{
			stats_capacity = 4;
			uint32_t old_stats[2] = { u_stats.emb_stats[0], u_stats.emb_stats[1] };

			u_stats.stats = new uint32_t[stats_capacity];
			u_stats.stats[0] = old_stats[0];
			u_stats.stats[1] = old_stats[1];
		}
		else
		{
			auto old_stats = u_stats.stats;

			switch (stats_capacity)
			{
			case 4: stats_capacity = 8;	break;
			case 8: stats_capacity = 12;	break;
			case 12: stats_capacity = 16;	break;
			case 16: stats_capacity = 20;	break;
			case 20: stats_capacity = 24;	break;
			case 24: stats_capacity = 32;	break;
			case 32: stats_capacity = 48;	break;
			case 48: stats_capacity = 64;	break;
			case 64: stats_capacity = 80;	break;
			case 80: stats_capacity = 96;	break;
			default:
				stats_capacity = (uint16_t)(stats_capacity * 1.4);
			}

			if (stats_capacity > compact_limit)
				stats_capacity = (uint16_t)(compact_limit + (uint16_t)1);

			u_stats.stats = new uint32_t[stats_capacity];
			copy_n(old_stats, stats_size, u_stats.stats);
			delete[] old_stats;
		}
	}

public:
	CAdjustableModelEmb() : total(0)
	{
		stats_capacity = 2;
		stats_size = 0;
	};

	~CAdjustableModelEmb()
	{
		if(stats_capacity != 2)
			delete[] u_stats.stats;
	};

	CAdjustableModelEmb(const CAdjustableModelEmb& c) = delete;
	CAdjustableModelEmb& operator=(const CAdjustableModelEmb&) = delete;

	// Warning: _init_stats are ignored in this model
	void Init(int* _init_stats)
	{
		if (stats_capacity != 2)
			delete[] u_stats.stats;

		total = NO_SYMBOLS;
		stats_size = 0;

		stats_capacity = 2;
	}

	void Init(const CAdjustableModelEmb& c)
	{
		total = c.total;

		if (stats_capacity != 2)
			delete[] u_stats.stats;

		stats_size = c.stats_size;
		stats_capacity = c.stats_capacity;

		if (stats_capacity == 2)
		{
			u_stats = c.u_stats;
		}
		else if (stats_capacity > 2)
		{
			u_stats.stats = new uint32_t[stats_capacity];
			copy_n(c.u_stats.stats, stats_size, u_stats.stats);
		}
		else
		{
			u_stats.stats = new uint32_t[NO_SYMBOLS];
			copy_n(c.u_stats.stats, NO_SYMBOLS, u_stats.stats);
		}
	}

	void GetFreq(int symbol, int& sym_freq, int& left_freq, int& totf)
	{
		left_freq = 0;

		if (stats_capacity == 2)
		{
			uint32_t cnt = 0;

			sym_freq = 1;

			for (uint32_t i = 0; i < stats_size; ++i)
			{
				auto x = u_stats.emb_stats[i];
				uint32_t s = x >> symbol_shift;
				uint32_t v = x & value_mask;

				if ((int)s < symbol)
				{
					left_freq += v;
					++cnt;
				}
				else if ((int)s == symbol)
				{
					sym_freq = v;
					break;
				}
				else
					break;
			}

			left_freq += symbol - cnt;
		}
		else if (stats_capacity)
		{
			uint32_t cnt = 0;

			sym_freq = 1;

			for (uint32_t i = 0; i < stats_size; ++i)
			{
				auto x = u_stats.stats[i];
				uint32_t s = x >> symbol_shift;
				uint32_t v = x & value_mask;

				if ((int)s < symbol)
				{
					left_freq += v;
					++cnt;
				}
				else if ((int)s == symbol)
				{
					sym_freq = v;
					break;
				}
				else
					break;
			}

			left_freq += symbol - cnt;
		}
		else
		{
			switch (symbol)
			{
			case 4: left_freq += u_stats.stats[3];
			case 3: left_freq += u_stats.stats[2];
			case 2: left_freq += u_stats.stats[1];
			case 1: left_freq += u_stats.stats[0];
			case 0: break;
			default:
				for (int i = 0; i < symbol; ++i)
					left_freq += u_stats.stats[i];
			}

			sym_freq = u_stats.stats[symbol];
		}

		totf = total;
	}

	void Update(int symbol)
	{
		if(stats_capacity == 2)
		{
			uint32_t i;

			for (i = 0; i < stats_size; ++i)
			{
				uint32_t s = u_stats.emb_stats[i] >> symbol_shift;

				if ((int)s == symbol)
					break;
				else if ((int)s > symbol)
				{
					if (stats_size == stats_capacity)
					{
						resize_stats();
						u_stats.stats[2] = u_stats.stats[1];
						if (i == 0)
							u_stats.stats[1] = u_stats.stats[0];
						
						u_stats.stats[i] = pack_sv(symbol, 1);
					}
					else 
					{
						u_stats.emb_stats[1] = u_stats.emb_stats[0];
						u_stats.emb_stats[0] = pack_sv(symbol, 1);
					}
					++stats_size;

					break;
				}
			}

			if (i == stats_size)
			{
				if (stats_size == stats_capacity)
				{
					resize_stats();
					u_stats.stats[stats_size] = pack_sv(symbol, 1);
				}
				else
					u_stats.emb_stats[stats_size] = pack_sv(symbol, 1);

				++stats_size;
			}

			if (stats_capacity == 2)
				u_stats.emb_stats[i] += ADDER;
			else
				u_stats.stats[i] += ADDER;
			total += ADDER;
		}
		else if (stats_capacity)
		{
			uint32_t i;
			bool expanded = false;

			for (i = 0; i < stats_size; ++i)
			{
				uint32_t s = u_stats.stats[i] >> symbol_shift;

				if ((int)s == symbol)
					break;
				else if ((int)s > symbol)
				{
					if (stats_size == stats_capacity)
						resize_stats();

					copy_backward(u_stats.stats + i, u_stats.stats + stats_size, u_stats.stats + stats_size + 1);
					u_stats.stats[i] = pack_sv(symbol, 1);
					++stats_size;

					expanded = true;
					break;
				}
			}

			if (i == stats_size)
			{
				if (stats_size == stats_capacity)
					resize_stats();

				u_stats.stats[stats_size] = pack_sv(symbol, 1);
				++stats_size;
				expanded = true;
			}

			u_stats.stats[i] += ADDER;
			total += ADDER;

			if (expanded && stats_size >= compact_limit)
			{
				auto old_stats = u_stats.stats;

				u_stats.stats = new uint32_t[NO_SYMBOLS];
				fill_n(u_stats.stats, NO_SYMBOLS, 1);
				for (uint32_t i = 0; i < stats_size; ++i)
				{
					auto x = old_stats[i];
					uint32_t s = x >> symbol_shift;
					uint32_t v = x & value_mask;

					u_stats.stats[s] = v;
				}

				stats_capacity = 0;

				delete[] old_stats;
			}
		}
		else
		{
			u_stats.stats[symbol] += ADDER;
			total += ADDER;
		}

		if (total >= MAX_TOTAL)
			rescale();
	}

	int GetSym(uint32_t left_freq)
	{
		uint32_t t = 0;

		if(stats_capacity == 2)
		{
			uint32_t cnt = 0;

			if(0 < stats_size)
			{
				auto x = u_stats.emb_stats[0];
				uint32_t s = x >> symbol_shift;
				uint32_t v = x & value_mask;

				t += v;

				if (t + (s - cnt) > left_freq)
				{
					if (t + (s - cnt) - v <= left_freq)
						return s;
					return s - (t + (s - cnt) - v - left_freq);
				}

				++cnt;
			}

			if(1 < stats_size)
			{
				auto x = u_stats.emb_stats[1];
				uint32_t s = x >> symbol_shift;
				uint32_t v = x & value_mask;

				t += v;

				if (t + (s - cnt) > left_freq)
				{
					if (t + (s - cnt) - v <= left_freq)
						return s;
					return s - (t + (s - cnt) - v - left_freq);
				}

				++cnt;
			}

			return NO_SYMBOLS - (total - left_freq);
		}
		else if (stats_capacity)
		{
			uint32_t cnt = 0;

			for (uint32_t i = 0; i < stats_size; ++i)
			{
				auto x = u_stats.stats[i];
				uint32_t s = x >> symbol_shift;
				uint32_t v = x & value_mask;

				t += v;

				if (t + (s - cnt) > left_freq)
				{
					if (t + (s - cnt) - v <= left_freq)
						return s;
					return s - (t + (s - cnt) - v - left_freq);
				}

				++cnt;
			}

			return NO_SYMBOLS - (total - left_freq);
		}
		else
		{
			for (uint32_t i = 0; i < NO_SYMBOLS; ++i)
			{
				t += u_stats.stats[i];
				if (t > left_freq)
					return i;
			}
		}

		return -1;
	}

	uint32_t GetTotal()
	{
		return total;
	}
};


// *******************************************************************************************
//
// *******************************************************************************************
template<typename T_MODEL, typename T_IO_STREAM, 
unsigned NO_SYMBOLS, unsigned MAX_LOG_COUNTER, unsigned ADDER> class CRangeCoderModel
{
	CRangeEncoder<T_IO_STREAM> *rce;
	CRangeDecoder<T_IO_STREAM> *rcd;

	T_MODEL model;

	const int rescale = 1u << MAX_LOG_COUNTER;
	bool compress;

public:
	CRangeCoderModel(CBasicRangeCoder<T_IO_STREAM> *rcb, int* _init, bool _compress) :
		compress(_compress)
	{
		model.Init(_init);

		if (compress)
		{
			rce = (CRangeEncoder<T_IO_STREAM>*) (rcb);
			rcd = nullptr;
		}
		else
		{
			rce = nullptr;
			rcd = (CRangeDecoder<T_IO_STREAM>*) (rcb);
		}
	}

	CRangeCoderModel(const CRangeCoderModel &c)
	{
		model.Init(c.model);
		rce = c.rce;
		rcd = c.rcd;

		compress = c.compress;
	}

	~CRangeCoderModel()
	{
	}

	void Encode(const int x)
	{
		int syfreq, ltfreq, totf;
		model.GetFreq(x, syfreq, ltfreq, totf);
		rce->EncodeFrequency(syfreq, ltfreq, totf);

		model.Update(x);
	}

	void Update(const int x)
	{
		model.Update(x);
	}

	int Decode()
	{
		int syfreq, ltfreq, totf;

		totf = model.GetTotal();
		ltfreq = (int) rcd->GetCumulativeFreq(totf);

		int x = model.GetSym(ltfreq);

		model.GetFreq(x, syfreq, ltfreq, totf);
		rcd->UpdateFrequency(syfreq, ltfreq, totf);
		model.Update(x);

		return x;
	}

	T_MODEL* GetSimpleModel()
	{
		return &model;
	}

	void Init(int *init)
	{
		model.Init(init);
	}
};

// *******************************************************************************************
//
// *******************************************************************************************
template<typename T_MODEL_FIXEXD_SIZE, typename T_IO_STREAM, unsigned N_SYMBOLS> class CRangeCoderModelFixedSize
{
	CRangeEncoder<T_IO_STREAM> *rce;
	CRangeDecoder<T_IO_STREAM> *rcd;

	T_MODEL_FIXEXD_SIZE model;

	bool compress;

public:
	CRangeCoderModelFixedSize(CBasicRangeCoder<T_IO_STREAM> *rcb, int _lg_totf, int _rescale, int* _init, uint32_t _adder, bool _compress) :
		compress(_compress)
	{
		model.Init(_init, _rescale, _adder);

		if (compress)
		{
			rce = (CRangeEncoder<T_IO_STREAM>*) (rcb);
			rcd = nullptr;
		}
		else
		{
			rce = nullptr;
			rcd = (CRangeDecoder<T_IO_STREAM>*) (rcb);
		}
	}

	CRangeCoderModelFixedSize(const CRangeCoderModelFixedSize &c)
	{
		model.Init(c.simple_model);
		rce = c.rce;
		rcd = c.rcd;

		compress = c.compress;
	}

	~CRangeCoderModelFixedSize()
	{
	}

	void Encode(const int x)
	{
		int syfreq, ltfreq, totf;

		model.GetFreq(x, syfreq, ltfreq, totf);
		rce->EncodeFrequency(syfreq, ltfreq, totf);

		model.Update(x);
	}

	int Decode()
	{
		int syfreq, ltfreq;
		int totf = model.GetTotal();
		ltfreq = (int) rcd->GetCumulativeFreq(totf);

		int x = model.GetSym(ltfreq);

		model.GetFreq(x, syfreq, ltfreq, totf);
		rcd->UpdateFrequency(syfreq, ltfreq, totf);
		model.Update(x);

		return x;
	}

	T_MODEL_FIXEXD_SIZE* GetSimpleModel()
	{
		return &model;
	}
};

// EOF
