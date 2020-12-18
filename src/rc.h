#pragma once
// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2020-12-18
// *******************************************************************************************

#include "defs.h"
#include "sub_rc.h"
#include <cmath>
#include <algorithm>
#include <numeric>

// *******************************************************************************************
//
// *******************************************************************************************
class CSimpleModel
{
	uint32_t n_symbols;
	uint32_t max_total;
	uint32_t *stats;
	uint32_t total;
	uint32_t adder;

	void rescale()
	{
		while (total >= max_total)
		{
			total = 0;
			for (uint32_t i = 0; i < n_symbols; ++i)
			{
				stats[i] = (stats[i] + 1) / 2;
				total += stats[i];
			}
		}
	}

public: 
	CSimpleModel(uint32_t _adder = 1) : n_symbols(0), max_total(0), stats(nullptr), total(0), adder(_adder)
	{
	};

	~CSimpleModel()
	{
		if (stats)
			delete[] stats;
	};

	CSimpleModel(const CSimpleModel &c) = delete;
	CSimpleModel& operator=(const CSimpleModel&) = delete;

	void Init(uint32_t _n_symbols, int *_init_stats, uint32_t _max_total, uint32_t _adder)
	{
		adder = _adder;

		if (stats)
		{
			if (n_symbols != _n_symbols)
			{
				delete[] stats;
				n_symbols = _n_symbols;
				stats = new uint32_t[n_symbols];
			}
		}
		else
		{
			n_symbols = _n_symbols;
			stats = new uint32_t[n_symbols];
		}

		max_total = _max_total;

		if (_init_stats)
			for (uint32_t i = 0; i < n_symbols; ++i)
				stats[i] = _init_stats[i];
		else
			fill_n(stats, n_symbols, 1);

		total = accumulate(stats, stats+n_symbols, 0u);
		rescale();
	}

	void Init(const CSimpleModel &c)
	{
		n_symbols = c.n_symbols;
		max_total = c.max_total;
		adder = c.adder;

		if (stats)
			delete[] stats;

		stats = new uint32_t[n_symbols];
		copy_n(c.stats, n_symbols, stats);
		total = accumulate(stats, stats + n_symbols, 0u);
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
		stats[symbol] += adder;
		total += adder;

		if (total >= max_total)
			rescale();
	}

	int GetSym(int left_freq)
	{
		int t = 0;

		for (uint32_t i = 0; i < n_symbols; ++i)
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
		for (uint32_t i = 0; i < n_symbols; ++i)
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
		for (uint32_t i = 0; i < n_symbols; ++i)
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
class CAdjustableModel
{
	const uint32_t value_mask = 0xffffffu;
	const uint32_t symbol_mask = 0xff000000u;
	const uint32_t symbol_shift = 24;
	const float compact_limit_frac = 0.25;

	uint32_t n_symbols;
	bool compact;
	uint32_t compact_limit;
	uint32_t max_total;
	vector<uint32_t> stats;
	uint32_t total;
	uint32_t adder;

	constexpr uint32_t pack_sv(uint32_t s, uint32_t v)
	{
		return (s << symbol_shift) + v;
	}

	void rescale()
	{
		while (total >= max_total)
		{
			total = 0;

			if (compact)
			{
				for (auto& x : stats)
				{
					uint32_t v = ((x & value_mask) + 1) / 2;
					x = (x & symbol_mask) + v;
					total += v;
				}

				total += n_symbols - (uint32_t) stats.size();
			}
			else
			{
				for (auto& x : stats)
				{
					x = (x + 1) / 2;
					total += x;
				}
			}
		}
	}

public:
	CAdjustableModel(uint32_t _adder = 1) : n_symbols(0), compact(true), max_total(0), total(0), adder(_adder)
	{
	};

	~CAdjustableModel()
	{
	};

	CAdjustableModel(const CAdjustableModel& c) = delete;
	CAdjustableModel& operator=(const CAdjustableModel&) = delete;

	// Warning: _init_stats are ignored in this model
	void Init(uint32_t _n_symbols, int* _init_stats, uint32_t _max_total, uint32_t _adder)
	{
		adder = _adder;
		n_symbols = _n_symbols;
		compact = true;
		max_total = _max_total;
		total = n_symbols;
		compact_limit = (uint32_t) ((double) n_symbols * compact_limit_frac);

		stats.clear();
		stats.shrink_to_fit();

	}

	void Init(const CAdjustableModel& c)
	{
		n_symbols = c.n_symbols;
		compact = c.compact;
		max_total = c.max_total;
		adder = c.adder;
		total = c.total;
		compact_limit = c.compact_limit;

		stats = c.stats;
	}

	void GetFreq(int symbol, int& sym_freq, int& left_freq, int& totf)
	{
		left_freq = 0;

		if (compact)
		{
			uint32_t cnt = 0;

			sym_freq = 1;

			for (auto x : stats)
			{
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
		if (compact)
		{
			uint32_t i;
			bool expanded = false;

			for (i = 0; i < stats.size(); ++i)
			{
				uint32_t s = stats[i] >> symbol_shift;

				if ((int) s == symbol)
					break;
				else if ((int) s > symbol)
				{
					if (stats.size() == stats.capacity())
					{
						if (stats.size() == 0)
							stats.reserve(1);
						else if (stats.size() < 3)
							stats.reserve(stats.size() * 2);
						else
							stats.reserve((size_t) ((double) stats.size() * 1.5));
					}

					stats.insert(stats.begin() + i, pack_sv(symbol, 1));
					expanded = true;
					break;
				}
			}

			if (i == stats.size())
			{
				stats.push_back(pack_sv(symbol, 1));
				expanded = true;
			}

			stats[i] += adder;
			total += adder;

			if (expanded && stats.size() > compact_limit)
			{
				vector<uint32_t> tmp(n_symbols, 1);

				for (auto x : stats)
				{
					uint32_t s = x >> symbol_shift;
					uint32_t v = x & value_mask;

					tmp[s] = v;
				}

				compact = false;
				stats = move(tmp);
			}
		}
		else
		{
			stats[symbol] += adder;
			total += adder;
		}

		if (total >= max_total)
			rescale();
	}

	int GetSym(uint32_t left_freq)
	{
		uint32_t t = 0;

		if (compact)
		{
			uint32_t cnt = 0;

			for (auto x : stats)
			{
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

			return n_symbols - (total - left_freq);
		}
		else
		{
			for (uint32_t i = 0; i < n_symbols; ++i)
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

	vector<uint32_t> GetStats()
	{
		return stats;
	}
};


// *******************************************************************************************
//
// *******************************************************************************************
template<typename T_MODEL, typename T_IO_STREAM> class CRangeCoderModel
{
	CRangeEncoder<T_IO_STREAM> *rce;
	CRangeDecoder<T_IO_STREAM> *rcd;

	T_MODEL model;

	int no_symbols;
	int lg_totf;
	int totf;
	int rescale;
	uint32_t adder;
	bool compress;

public:
	CRangeCoderModel(CBasicRangeCoder<T_IO_STREAM> *rcb, int _no_symbols, int _lg_totf, int _rescale, int* _init, uint32_t _adder, bool _compress) :
		no_symbols(_no_symbols), lg_totf(_lg_totf), totf(1 << _lg_totf), rescale(_rescale), adder(_adder), compress(_compress)
	{
		model.Init(no_symbols, _init, rescale, adder);

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

		no_symbols = c.no_symbols;
		lg_totf = c.lg_totf;
		totf = c.totf;
		rescale = c.rescale;
		compress = c.compress;
		adder = c.adder;
	}

	~CRangeCoderModel()
	{
	}

	void Encode(const int x)
	{
		int syfreq, ltfreq;
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
		int syfreq, ltfreq;

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
		model.Init(no_symbols, init, rescale, adder);
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
	double est_tot_len;

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

		est_tot_len = 0;
	}

	CRangeCoderModelFixedSize(const CRangeCoderModelFixedSize &c)
	{
		model.Init(c.simple_model);
		rce = c.rce;
		rcd = c.rcd;

		compress = c.compress;

		est_tot_len = 0;
	}

	~CRangeCoderModelFixedSize()
	{
	}

	void Encode(const int x)
	{
		int syfreq, ltfreq, totf;

		model.GetFreq(x, syfreq, ltfreq, totf);
		rce->EncodeFrequency(syfreq, ltfreq, totf);

		est_tot_len += rce->EstimateCodeLen(syfreq, totf);

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

	double GetEstTotLen()
	{
		return est_tot_len;
	}
};

// EOF
