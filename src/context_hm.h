#pragma once
// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include <mmintrin.h>
#include <cstdint>
#include <xmmintrin.h>
#include <iostream> 
#include <cstddef>

#include "defs.h"
#include "rc.h"
#include "io.h"

//#define ENABLE_CTX_HM_COUNTER

// ************************************************************************************
template<typename MODEL> class CContextHM {
public:
	typedef struct {
		context_t ctx;
		MODEL *rcm;
#ifdef ENABLE_CTX_HM_COUNTER
		size_t counter;
#endif
	} item_t;

	typedef context_t key_type;
	typedef MODEL* value_type;
	typedef size_t aux_type;

private:
	double max_fill_factor;

	size_t size;
	item_t *data;
	size_t allocated;
	size_t size_when_restruct;
	size_t allocated_mask;

	size_t ht_memory;
	size_t ht_total;
	size_t ht_match;

	void restruct(void)
	{
		item_t *old_data = data;
		size_t old_allocated = allocated;

		allocated *= 2;
		size = 0;

		allocated_mask = allocated - 1ull;
		size_when_restruct = (size_t)((double) allocated * max_fill_factor);

		data = new item_t[allocated];
		for (size_t i = 0; i < allocated; ++i)
			data[i].rcm = nullptr;

		ht_memory += allocated * sizeof(item_t);

		for (size_t i = 0; i < old_allocated; ++i)
			if (old_data[i].rcm != nullptr)
#ifdef ENABLE_CTX_HM_COUNTER
				insert(old_data[i].ctx, old_data[i].rcm, old_data[i].counter);
#else
				insert(old_data[i].ctx, old_data[i].rcm);
#endif

		delete[] old_data;
		ht_memory -= old_allocated * sizeof(item_t);
	}

	// Based on murmur64
	size_t hash(context_t ctx)
	{
		auto h = ctx;

		h ^= h >> 33;
		h *= 0xff51afd7ed558ccdL;
		h ^= h >> 33;
		h *= 0xc4ceb9fe1a85ec53L;
		h ^= h >> 33;
		
		return h & allocated_mask;
	}

public:
	CContextHM()
	{
		ht_memory = 0;
		ht_total = 0;
		ht_match = 0;

		allocated = 1u << 5;
		allocated_mask = allocated - 1;

		size = 0;
		data = new item_t[allocated];
		for (size_t i = 0; i < allocated; ++i)
			data[i].rcm = nullptr;

		max_fill_factor = 0.6;

		ht_memory += allocated * sizeof(item_t);

		size_when_restruct = (size_t)((double) allocated * max_fill_factor);
	}

	~CContextHM()
	{
		if (data == nullptr)
			return;

		for (size_t i = 0; i < allocated; ++i)
			if (data[i].rcm)
				delete data[i].rcm;
		delete[] data;
	}

	size_t get_bytes() const {
		return ht_memory;
	}

	void debug_list(vector<CContextHM<MODEL>::item_t> &v_ctx)
	{
		v_ctx.clear();

		for (size_t i = 0; i < allocated; ++i)
			if (data[i].rcm)
				v_ctx.push_back(data[i]);

		sort(v_ctx.begin(), v_ctx.end(), [](auto &x, auto &y) {return x.counter > y.counter; });
	}

#ifdef ENABLE_CTX_HM_COUNTER
	bool insert(const context_t ctx, MODEL *rcm, size_t counter = 0)
#else
	bool insert(const context_t ctx, MODEL *rcm)
#endif
	{
		if (size >= size_when_restruct)
			restruct();

		size_t h = hash(ctx);

		if (data[h].rcm != nullptr)
		{
			do
			{
				h = (h + 1) & allocated_mask;
			} while (data[h].rcm != nullptr);
		}

		++size;

		data[h].ctx = ctx;
		data[h].rcm = rcm;
#ifdef ENABLE_CTX_HM_COUNTER
		data[h].counter = counter;
#endif

		return true;
	}

	MODEL* find(const context_t ctx)
	{
		size_t h = hash(ctx);

		if (data[h].rcm == nullptr)
			return nullptr;

		if (data[h].ctx == ctx)
			return data[h].rcm;

		h = (h + 1) & allocated_mask;

		while (data[h].rcm != nullptr)
		{
			if (data[h].ctx == ctx)
				return data[h].rcm;
			h = (h + 1) & allocated_mask;
		}

		return nullptr;
	}

#ifdef ENABLE_CTX_HM_COUNTER
	MODEL* find_ext(const context_t ctx, size_t *&p_counter)
	{
		size_t h = hash(ctx);

		if (data[h].rcm == nullptr)
			return nullptr;
		
		if (data[h].ctx == ctx)
		{
			p_counter = &data[h].counter;
			return data[h].rcm;
		}

		h = (h + 1) & allocated_mask;

		while (data[h].rcm != nullptr)
		{
			if (data[h].ctx == ctx)
			{
				p_counter = &data[h].counter;
				return data[h].rcm;
			}
			h = (h + 1) & allocated_mask;
		}

		return nullptr;
	}
#endif

	void prefetch(const context_t ctx)
	{
		size_t h = hash(ctx);

#ifdef _WIN32
		_mm_prefetch((const char*)(data + h), _MM_HINT_T0);
#else
		__builtin_prefetch(data + h);
#endif
	}

	size_t get_size(void) const
	{
		return size;
	}
}; 

// EOF
