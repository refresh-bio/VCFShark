#pragma once
// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include <algorithm>
#include <vector>
#include <cstdint>
#include <cstring>

#include <iostream>

using namespace std;

#ifndef _WIN32
#define my_fseek	fseek
#define my_ftell	ftell
#else
#define my_fseek	_fseeki64
#define my_ftell	_ftelli64
#endif

// *******************************************************************************************
// Buffered input file
class CInFile
{
	const size_t BUFFER_SIZE = 128 << 20;

	FILE *f;
	uint8_t *buffer;
	size_t buffer_pos;
	size_t buffer_filled;

	size_t file_size;
	size_t before_buffer_bytes;

public:
	CInFile() : f(nullptr), buffer(nullptr)
	{};

	~CInFile()
	{
		if (f)
			fclose(f);
		if (buffer)
			delete[] buffer;
	}

	bool Open(string file_name)
	{
		if (f)
			return false;

		f = fopen(file_name.c_str(), "rb");
		if (!f)
			return false;

		my_fseek(f, 0, SEEK_END);
		file_size = my_ftell(f);
		my_fseek(f, 0, SEEK_SET);
		before_buffer_bytes = 0;

		buffer = new uint8_t[BUFFER_SIZE];
		buffer_pos = 0;
		buffer_filled = 0;

#ifdef LOG_INFO
		cout << "Opening file of size: " << file_size << endl;
#endif

		return true;
	}

	bool Close()
	{
		if (f)
		{
			fclose(f);
			f = nullptr;
		}
		if (buffer)
		{
			delete[] buffer;
			buffer = nullptr;
		}

		return true;
	}

	int Get()
	{
		if (buffer_pos < buffer_filled)
			return buffer[buffer_pos++];

		if (feof(f))
			return EOF;

		before_buffer_bytes += buffer_filled;

		buffer_filled = fread(buffer, 1, BUFFER_SIZE, f);
		if (buffer_filled == 0)
			return EOF;

		buffer_pos = 0;
		return buffer[buffer_pos++];
	}

	uint8_t GetByte()
	{
		return (uint8_t)Get();
	}

	uint64_t ReadUInt(int no_bytes)
	{
		uint64_t x = 0;
		uint64_t shift = 0;

		for (int i = 0; i < no_bytes; ++i)
		{
			uint64_t c = Get();
			x += c << shift;
			shift += 8;
		}

		return x;
	}

	// !!! To wygl¹da na b³¹d - brak sprawdzania 
	void Read(uint8_t *ptr, uint64_t size)
	{
		if (before_buffer_bytes + buffer_pos + size > file_size)
			size = file_size - (before_buffer_bytes + buffer_pos);

		uint64_t to_read = size;

		while (buffer_pos + to_read > BUFFER_SIZE)
		{
			memcpy(ptr, buffer + buffer_pos, BUFFER_SIZE - buffer_pos);
			ptr += BUFFER_SIZE - buffer_pos;
			to_read -= BUFFER_SIZE - buffer_pos;

			before_buffer_bytes += buffer_filled;
			buffer_filled = fread(buffer, 1, BUFFER_SIZE, f);
			buffer_pos = 0;
		}

		memcpy(ptr, buffer + buffer_pos, to_read);
		buffer_pos += to_read;
	}

	bool Eof()
	{
		return before_buffer_bytes + buffer_pos >= file_size;
	}

	size_t FileSize()
	{
		if (f)
			return file_size;
		else
			return 0;
	}

	size_t GetPos()
	{
		return before_buffer_bytes + buffer_pos;
	}
};

// *******************************************************************************************
// Buffered output file
class COutFile
{
	const size_t BUFFER_SIZE = 8 << 18; //8 << 20;

	FILE *f;
	uint8_t *buffer;
	size_t buffer_pos;
	bool success;

public:
	COutFile() : f(nullptr), buffer(nullptr)
	{};

	~COutFile()
	{
		if (f)
			Close();
		if (buffer)
			delete[] buffer;
	}

	bool Open(string file_name)
	{
		if (f)
			return false;

		f = fopen(file_name.c_str(), "wb");
		if (!f)
			return false;

		buffer = new uint8_t[BUFFER_SIZE];
		buffer_pos = 0;
		success = true;

		return true;
	}

	bool Close()
	{
		if (!f)
			return true;

		if (buffer_pos)
			success &= fwrite(buffer, 1, buffer_pos, f) == buffer_pos;

		if (f)
		{
			fclose(f);
			f = nullptr;
		}
		if (buffer)
		{
			delete[] buffer;
			buffer = nullptr;
		}

		return success;
	}

	void PutByte(uint8_t c)
	{
		if (buffer_pos == BUFFER_SIZE)
		{
			success &= fwrite(buffer, 1, BUFFER_SIZE, f) == BUFFER_SIZE;
			buffer_pos = 0;
		}

		buffer[buffer_pos++] = c;
	}

	void Put(char c)
	{
		if (buffer_pos == BUFFER_SIZE)
		{
			success &= fwrite(buffer, 1, BUFFER_SIZE, f) == BUFFER_SIZE;
			buffer_pos = 0;
		}

		buffer[buffer_pos++] = c;
	}

	void Write(const uint8_t *p, size_t n)
	{
		uint8_t *q = (uint8_t *)p;

		while (buffer_pos + n > BUFFER_SIZE)
		{
			size_t small_n = BUFFER_SIZE - buffer_pos;
			memcpy(buffer + buffer_pos, q, small_n);
			success &= fwrite(buffer, 1, BUFFER_SIZE, f) == BUFFER_SIZE;

			buffer_pos = 0;
			n -= small_n;
			q += small_n;
		}

		memcpy(buffer + buffer_pos, q, n);
		buffer_pos += n;
	}

	void WriteUInt(uint64_t x, int no_bytes)
	{
		for (int i = 0; i < no_bytes; ++i)
		{
			PutByte(x & 0xffu);
			x >>= 8;
		}
	}

	void Write(string &s)
	{
		Write((uint8_t*)s.c_str(), s.size());
	}

	void Write(string &s, size_t start_pos, size_t len)
	{
		Write((uint8_t*)s.c_str() + start_pos, len);
	}
};


// *******************************************************************************************
// Class for storage of range coder compressed data
class CVectorIOStream
{
	vector<uint8_t> &v;
	size_t read_pos;

public:
	CVectorIOStream(vector<uint8_t> &_v) : v(_v), read_pos(0)
	{}

	void RestartRead()
	{
		read_pos = 0;
	}

	bool Eof()
	{
		return read_pos >= v.size();
	}

	uint8_t GetByte()
	{
		return v[read_pos++];
	}

	void PutByte(uint8_t x)
	{
		v.push_back(x);
	}

	size_t Size()
	{
		return v.size();
	}

	uint8_t* Data()
	{
		return v.data();
	}
};

// EOF
