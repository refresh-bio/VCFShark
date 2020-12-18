#pragma once
// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2020-12-18
// *******************************************************************************************

#ifndef H_RANGECODER
#define H_RANGECODER

#include <iostream>
#include <cmath>
#include <cstdint>

#include "defs.h"
#include <assert.h>


// *******************************************************************************************
//
// *******************************************************************************************
template<typename T_IO_STREAM> class CBasicRangeCoder
{
};

// *******************************************************************************************
//
// *******************************************************************************************
template<typename T_IO_STREAM> class CRangeEncoder : public CBasicRangeCoder<T_IO_STREAM>
{
public:
	typedef uint64_t Code;
	typedef uint64_t Freq;

	static const Freq TopValue = 0x00ffffffffffffULL;
	static const Code Mask64 = 0xff00000000000000ULL;

	T_IO_STREAM &io_stream;

	Code	low;
	Freq	range;

	CRangeEncoder(T_IO_STREAM& _io_stream) : io_stream(_io_stream), low(0), range(0)
	{}

	void Start()
	{
		low = 0;
		range = Mask64;
	}

	float EstimateCodeLen(Freq s, Freq t)
	{
		return -log2f((float)s / t);
	}

	void EncodeFrequency(Freq symFreq_, Freq cumFreq_, Freq totalFreqSum_)
	{
		assert(range > totalFreqSum_);
		range /= totalFreqSum_;
		low += range * cumFreq_;
		range *= symFreq_;

		while (range <= TopValue)
		{
			assert(range != 0);
			if ((low ^ (low + range)) & Mask64)
			{
				Freq r = (Freq)low;
				range = (r | TopValue) - r;
			}
			io_stream.PutByte((uint8_t) (low >> 56));
			low <<= 8, range <<= 8;
		}
	}

	void End()
	{
		for (int i = 0; i < 8; i++)
		{
			io_stream.PutByte((uint8_t) (low >> 56));
			low <<= 8;
		}
	}
};

// *******************************************************************************************
//
// *******************************************************************************************
template<typename T_IO_STREAM> class CRangeDecoder : public CBasicRangeCoder<T_IO_STREAM>
{
public:
	typedef uint64_t Code;
	typedef uint64_t Freq;

	static const Freq TopValue = 0x00ffffffffffffULL;
	static const Code Mask64 = 0xff00000000000000ULL;

	T_IO_STREAM &io_stream;

	Code	low;
	Freq	range;

	CRangeDecoder(T_IO_STREAM &_io_stream) : io_stream(_io_stream), low(0), range(0)
	{
		buffer = 0;
	}

	void Start()
	{
		if (io_stream.Size() < 8)
			return;

		buffer = 0;
		for (uint32_t i = 1; i <= 8; ++i)
		{
			buffer |= (Code)io_stream.GetByte() << (64 - i * 8);
		}

		low = 0;
		range = Mask64;
	}

	Freq GetCumulativeFreq(Freq totalFreq_)
	{
		assert(totalFreq_ != 0);
		return (Freq) (buffer / (range /= totalFreq_));
	}

	void UpdateFrequency(Freq symFreq_, Freq lowEnd_, Freq /*totalFreq_*/)
	{
		Freq r = lowEnd_*range;
		buffer -= r;
		low += r;
		range *= symFreq_;

		while (range <= TopValue)
		{
			if ((low ^ (low + range)) & Mask64)
			{
				Freq r = (Freq)low;
				range = (r | TopValue) - r;
			}

			buffer = (buffer << 8) + io_stream.GetByte();
			low <<= 8, range <<= 8;
		}
	}

	void End()
	{}

private:
	Code buffer;
};

#endif 

// EOF
