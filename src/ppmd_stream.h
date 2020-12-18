#pragma once
#include "common.h"
#include "coro3b_fake.h"
#include "libpmd.h"

#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
using namespace std;

// *******************************************************************************************
//
// *******************************************************************************************
class CPPMDStream
{
	ALIGN(4096) pmd_codec codec;

	size_t kernel_buf_size;
	size_t kernel_allocated_size;

	vector<uint8_t> kernel_input_buf;
	vector<uint8_t> kernel_output_buf;

	vector<uint8_t> raw_buf;
	vector<uint8_t> compressed_buf;

	size_t raw_pos;
	size_t raw_size;
	size_t uncompressed_size;

	bool is_encoding;

	bool init(bool _is_encoding, size_t order, size_t memory, size_t restore, size_t file_size);

	size_t size_with_overhead(size_t x)
	{
		return (size_t)(x * 1.2 + 64);
	}

	void kernel_resize(size_t _kernel_buf_size)
	{
		kernel_buf_size = _kernel_buf_size;
		kernel_allocated_size = size_with_overhead(kernel_buf_size);

		kernel_input_buf.reserve(kernel_allocated_size);
		kernel_output_buf.reserve(kernel_allocated_size);

		raw_size = 0;
	}

	void write_bytes(uint8_t* p, size_t n, bool check_space = true);
	void read_bytes(uint8_t* p, size_t n, bool check_space = true);
	char read_char();

	void compress();
	void decompress();

public:
	CPPMDStream();
	~CPPMDStream();

	bool InitCompress(size_t _kernel_buf_size = 32, unsigned int order = 12, unsigned int memory = 256, unsigned int restore = 1);
	bool InitDecompress(size_t _kernel_buf_size = 32, unsigned int order = 12, unsigned int memory = 256, unsigned int restore = 1);

	void WriteUint64(uint64_t x);
	void WriteUint32(uint32_t x);
	void WriteInt64(uint64_t x);
	void WriteInt32(uint32_t x);
	void WriteString(string& s);
	void WriteArray(uint8_t *p, uint32_t n);
	void WriteStream(vector<uint8_t>& stream);

	uint64_t ReadUint64();
	uint32_t ReadUint32();
	int64_t ReadInt64();
	int32_t ReadInt32();
	string ReadString();
	void ReadArray(uint8_t* &p, uint32_t &n);
	void ReadStream(vector<uint8_t>& stream);

	size_t FinishCompress();
	bool CompleteDecompress();

	void SetCompressedVector(vector<uint8_t>& vec, size_t raw_size);
	void GetCompressedVector(vector<uint8_t>& vec, size_t &raw_size);
	bool NeedCompressedVector();
	bool FullCompressedVector();
};



// *******************************************************************************************
CPPMDStream::CPPMDStream()
{
	//raw_buf.reserve(kernel_allocated_size);
	//compressed_buf.reserve(kernel_allocated_size);
}

// *******************************************************************************************
CPPMDStream::~CPPMDStream()
{
}

// *******************************************************************************************
bool CPPMDStream::InitCompress(size_t _kernel_buf_size, unsigned int order, unsigned int memory, unsigned int restore)
{
	init(true, order, memory, restore, ~0ull);

	kernel_resize(_kernel_buf_size);

	return true;
}

// *******************************************************************************************
bool CPPMDStream::InitDecompress(size_t _kernel_buf_size, unsigned int order, unsigned int memory, unsigned int restore)
{
	init(false, order, memory, restore, ~0ull);

	kernel_resize(_kernel_buf_size);

	return true;
}

// *******************************************************************************************
size_t CPPMDStream::FinishCompress()
{
	compress();

	return raw_size;
}

// *******************************************************************************************
bool CPPMDStream::CompleteDecompress()
{
	return false;		// !! To remove
}

// *******************************************************************************************
void CPPMDStream::SetCompressedVector(vector<uint8_t>& vec, size_t raw_size)
{
	compressed_buf = move(vec);
	uncompressed_size = raw_size;

	kernel_resize(raw_size);

	decompress();
}

// *******************************************************************************************
void CPPMDStream::GetCompressedVector(vector<uint8_t>& vec, size_t& raw_size)
{
	vec = move(compressed_buf);
	raw_size = uncompressed_size;

	compressed_buf.clear();
	uncompressed_size = 0;
}

// *******************************************************************************************
bool CPPMDStream::NeedCompressedVector()
{
	return raw_buf.empty();
}

// *******************************************************************************************
bool CPPMDStream::FullCompressedVector()
{
	return !compressed_buf.empty();
}

// *******************************************************************************************
bool CPPMDStream::init(bool _is_encoding, size_t order, size_t memory, size_t restore, size_t file_size)
{
	uint32_t pmd_args[4];

	pmd_args[0] = order;
	pmd_args[1] = memory;
	pmd_args[2] = restore;
	pmd_args[3] = file_size;
	is_encoding = _is_encoding;

	cout << pmd_args[0] << " ";
	cout << pmd_args[1] << " ";
	cout << pmd_args[2] << " ";
	cout << pmd_args[3] << " ";
	cout << is_encoding << endl;


	auto x = codec.Init(!is_encoding, pmd_args);

	cout << "x: " << x << endl;

	if (x)
	{
		cout << "x: " << x << endl;
		cerr << "Cannot initialize PPMD\n";
		exit(1);
	}
}

// *******************************************************************************************
void CPPMDStream::compress()
{
	uncompressed_size = raw_buf.size();
	kernel_input_buf = move(raw_buf);
	raw_buf.clear();

	size_t kernel_output_size = size_with_overhead(uncompressed_size);
	kernel_output_buf.resize(kernel_output_size);

	codec.rc_Init();

	for (int i = 0; i < uncompressed_size; ++i)
		codec.EncodeByte(kernel_input_buf[i]);

	codec.rc_Quit();

	codec.get_output(kernel_output_buf);

	auto compressed_size = codec.getoutsize();
	compressed_buf.assign(kernel_output_buf.begin(), kernel_output_buf.begin() + compressed_size);
}

// *******************************************************************************************
void CPPMDStream::decompress()
{
	kernel_input_buf = move(compressed_buf);
	compressed_buf.clear();

	raw_buf.resize(uncompressed_size);


	codec.add_input(kernel_input_buf);

	codec.rc_Init();
	for (int i = 0; i < uncompressed_size; ++i)
		raw_buf[i] = codec.DecodeByte();

	codec.rc_Quit();

	compressed_buf.clear();

	raw_pos = 0;
}

// *******************************************************************************************
void CPPMDStream::write_bytes(uint8_t* p, size_t n, bool check_space)
{
	raw_buf.insert(raw_buf.end(), p, p + n);
	raw_size += n; 

	if (check_space && raw_buf.size() >= kernel_buf_size)
		compress();
}

// *******************************************************************************************
void CPPMDStream::read_bytes(uint8_t* p, size_t n, bool check_space)
{
	if (raw_buf.empty())
		decompress();

	copy_n(raw_buf.begin() + raw_pos, n, p);
	raw_pos += n;

	if (check_space && raw_pos == raw_buf.size())
	{
		raw_buf.clear();
		raw_pos = 0;
	}

	raw_size += n;
}

// *******************************************************************************************
char CPPMDStream::read_char()
{
	if (raw_buf.empty())
		decompress();

	char c = (char)raw_buf[raw_pos++];

	if (raw_pos == raw_buf.size())
	{
		raw_buf.clear();
		raw_pos = 0;
	}

	++raw_size;

	return c;
}

// *******************************************************************************************
void CPPMDStream::WriteUint64(uint64_t x)
{
	write_bytes((uint8_t*)&x, 8);
}

// *******************************************************************************************
void CPPMDStream::WriteUint32(uint32_t x)
{
	write_bytes((uint8_t*)&x, 4);
}

// *******************************************************************************************
void CPPMDStream::WriteInt64(uint64_t x)
{
	write_bytes((uint8_t*)&x, 8);
}

// *******************************************************************************************
void CPPMDStream::WriteInt32(uint32_t x)
{
	write_bytes((uint8_t*)&x, 4);
}

// *******************************************************************************************
void CPPMDStream::WriteString(string& s)
{
	write_bytes((uint8_t*)s.data(), s.size() + 1);
}

// *******************************************************************************************
void CPPMDStream::WriteArray(uint8_t* p, uint32_t n)
{
	write_bytes((uint8_t*)&n, 4, false);
	write_bytes(p, n);
}

// *******************************************************************************************
void CPPMDStream::WriteStream(vector<uint8_t>& stream)
{
	kernel_resize(stream.size());

	uint32_t stream_size = stream.size();

	write_bytes(stream.data(), stream_size, false);
	compress();
}

// *******************************************************************************************
uint64_t CPPMDStream::ReadUint64()
{
	uint64_t x;

	read_bytes((uint8_t*)&x, 8);

	return x;
}

// *******************************************************************************************
uint32_t CPPMDStream::ReadUint32()
{
	uint32_t x;

	read_bytes((uint8_t*)&x, 4);

	return x;
}

// *******************************************************************************************
int64_t CPPMDStream::ReadInt64()
{
	int64_t x;

	read_bytes((uint8_t*)&x, 8);

	return x;
}

// *******************************************************************************************
int32_t CPPMDStream::ReadInt32()
{
	int32_t x;

	read_bytes((uint8_t*)&x, 4);

	return x;
}

// *******************************************************************************************
string CPPMDStream::ReadString()
{
	string s;

	uint8_t c;

	while (true)
	{
		c = read_char();

		if (c != 0)
			s.push_back(c);
		else
			break;
	}

	return s;
}

// *******************************************************************************************
void CPPMDStream::ReadArray(uint8_t*& p, uint32_t& n)
{
	read_bytes((uint8_t*)&n, 4, false);

	if (p)
		delete[] p;
	p = new uint8_t[n];
	read_bytes(p, n);
}

// *******************************************************************************************
void CPPMDStream::ReadStream(vector<uint8_t>& stream)
{
	stream.resize(uncompressed_size);
	read_bytes(stream.data(), uncompressed_size, false);
}

// EOF
