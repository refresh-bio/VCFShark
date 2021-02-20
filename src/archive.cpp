// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include "archive.h"

#include <iostream>
#include <cstdio>
#include <utility>

#ifndef _WIN32
#define my_fseek	fseek
#define my_ftell	ftell
#else
#define my_fseek	_fseeki64
#define my_ftell	_ftelli64
#endif

using namespace std;

// ******************************************************************************
CArchive::CArchive(bool _input_mode)
{
	f = nullptr;
	input_mode = _input_mode;
}

// ******************************************************************************
CArchive::~CArchive()
{
	if (f)
		Close();
}

// ******************************************************************************
bool CArchive::Open(string _file_name)
{
	lock_guard<mutex> lck(mtx);

	if (f)
		fclose(f);

	m_streams.clear();
	file_name = _file_name;

	f = fopen(file_name.c_str(), input_mode ? "rb" : "wb");

	setvbuf(f, nullptr, _IOFBF, 64 << 20);

	if (!f)
		return false;

	if (input_mode)
		deserialize();

	f_offset = 0;

	return true;
}

// ******************************************************************************
bool CArchive::Close()
{
	lock_guard<mutex> lck(mtx);

	if (!f)
		return false;

	if (input_mode)
	{
		fclose(f);
		f = nullptr;
	}
	else
	{
		serialize();
		fclose(f);
		f = nullptr;
	}

	return true;
}

// ******************************************************************************
size_t CArchive::write_fixed(size_t x, FILE* file)
{
	fwrite(&x, 1, 8, file);

	return 8;
}

// ******************************************************************************
size_t CArchive::write(size_t x, FILE* file)
{
	int no_bytes = 0;

	for (size_t tmp = x; tmp; tmp >>= 8)
		++no_bytes;
	
	putc(no_bytes, file);

	for (int i = no_bytes; i; --i)
		putc((x >> ((i - 1) * 8)) & 0xff, file);

	return no_bytes + 1;
}

// ******************************************************************************
size_t CArchive::write(string s, FILE* file)
{
	fwrite(s.c_str(), 1, s.size(), file);
	putc(0, file);

	return s.size() + 1;
}

// ******************************************************************************
size_t CArchive::read_fixed(size_t& x, FILE* file)
{
	return fread(&x, 8, 1, file) * 8;
}

// ******************************************************************************
size_t CArchive::read(size_t& x, FILE* file)
{
	int no_bytes = getc(file);

	x = 0;

	for (int i = 0; i < no_bytes; ++i)
	{
		x <<= 8;
		x += (size_t)getc(file);
	}

	return no_bytes + 1;
}

// ******************************************************************************
size_t CArchive::read(string& s, FILE* file)
{
	s.clear();

	while (true)
	{
		int c = getc(file);
		if (c == EOF)
			return 0;

		if (c == 0)
			return s.size() + 1;

		s.push_back((char)c);
	}

	return 0;
}

// ******************************************************************************
bool CArchive::serialize()
{
	size_t footer_size = 0;

	// Store stream part offsets
	footer_size += write(m_streams.size(), f);

	for (auto& stream : m_streams)
	{
		size_t str_size = 0;

		footer_size += write(stream.second.stream_name, f);
		footer_size += write(stream.second.parts.size(), f);
		footer_size += write(stream.second.raw_size, f);

		for (auto& part : stream.second.parts)
		{
			footer_size += write(part.offset, f);
			footer_size += write(part.size, f);

			str_size += part.size;
		}

#ifdef LOG_INFO
		cerr << stream.first << ": " << stream.second.stream_name << "  raw size: " << stream.second.raw_size << "   packed size: " << str_size << endl;
#endif
	}

	write_fixed(footer_size, f);

	return true;
}

// ******************************************************************************
bool CArchive::deserialize()
{
	size_t footer_size;

	my_fseek(f, -8, SEEK_END);
	read_fixed(footer_size, f);

	my_fseek(f, -(long)(8 + footer_size), SEEK_END);

	// Load stream part offsets
	size_t n_streams;
	read(n_streams, f);

	for (size_t i = 0; i < n_streams; ++i)
	{
		m_streams[(int) i] = stream_t();
		auto& stream_second = m_streams[(int) i];

		read(stream_second.stream_name, f);
		read(stream_second.cur_id, f);
		read(stream_second.raw_size, f);

		stream_second.parts.resize(stream_second.cur_id);
		for (size_t j = 0; j < stream_second.cur_id; ++j)
		{
			read(stream_second.parts[j].offset, f);
			read(stream_second.parts[j].size, f);
		}

		stream_second.cur_id = 0;
	}
	
	my_fseek(f, 0, SEEK_SET);

	return true;
}

// ******************************************************************************
int CArchive::RegisterStream(string stream_name)
{
	lock_guard<mutex> lck(mtx);

	int id = (int) m_streams.size();

	m_streams[id] = stream_t();
	m_streams[id].cur_id = 0;
	m_streams[id].stream_name = stream_name;

	return id;
}

// ******************************************************************************
int CArchive::GetStreamId(string stream_name)
{
	lock_guard<mutex> lck(mtx);

	for (auto& x : m_streams)
		if (x.second.stream_name == stream_name)
			return x.first;

	return -1;
}

// ******************************************************************************
bool CArchive::AddPart(int stream_id, vector<uint8_t> &v_data, size_t metadata)
{
	lock_guard<mutex> lck(mtx);

	m_streams[stream_id].parts.emplace_back(f_offset, v_data.size());
	m_streams[stream_id].signatures.emplace_back(signature(v_data) ^ metadata);

	f_offset += write(metadata, f);

	if(v_data.size())
		fwrite(v_data.data(), 1, v_data.size(), f);

	f_offset += v_data.size();

	return true;
}

// ******************************************************************************
int CArchive::AddPartPrepare(int stream_id)
{
	lock_guard<mutex> lck(mtx);
	
	m_streams[stream_id].parts.emplace_back(0, 0);
	m_streams[stream_id].signatures.emplace_back(0);

	return (int) m_streams[stream_id].parts.size() - 1;
}

// ******************************************************************************
bool CArchive::AddPartComplete(int stream_id, int part_id, vector<uint8_t>& v_data, size_t metadata)
{
	lock_guard<mutex> lck(mtx);
	
	m_streams[stream_id].parts[part_id] = part_t(f_offset, v_data.size());
	m_streams[stream_id].signatures[part_id] = signature(v_data) ^ metadata;

	f_offset += write(metadata, f);
	if (v_data.size())
		fwrite(v_data.data(), 1, v_data.size(), f);

	f_offset += v_data.size();

	return true;
}

// ******************************************************************************
void CArchive::SetRawSize(int stream_id, size_t raw_size)
{
	lock_guard<mutex> lck(mtx);
	
	m_streams[stream_id].raw_size = raw_size;
}

// ******************************************************************************
size_t CArchive::GetRawSize(int stream_id)
{
	lock_guard<mutex> lck(mtx);
	
	return m_streams[stream_id].raw_size;
}

// ******************************************************************************
size_t CArchive::GetCompressedSize(int stream_id)
{
	lock_guard<mutex> lck(mtx);

	size_t size = 0;

	auto p = m_streams.find(stream_id);

	if (p == m_streams.end())
		return 0;

	for (auto q : p->second.parts)
		size += q.size;

	return size;
}

// ******************************************************************************
bool CArchive::GetPart(int stream_id, vector<uint8_t> &v_data, size_t &metadata)
{
	lock_guard<mutex> lck(mtx);
	
	auto& p = m_streams[stream_id];

	if (p.cur_id >= p.parts.size())
		return false;

	v_data.resize(p.parts[p.cur_id].size);

	my_fseek(f, p.parts[p.cur_id].offset, SEEK_SET);

	if(p.parts[p.cur_id].size != 0)
		read(metadata, f);
	else
	{
		metadata = 0;
		p.cur_id++;
		return true;
	}

	auto r = fread(v_data.data(), 1, p.parts[p.cur_id].size, f);

	p.cur_id++;

	if (r != p.parts[p.cur_id-1].size)
		return false;

	return r == p.parts[p.cur_id-1].size;
}

// ******************************************************************************
bool CArchive::ResetStreamPartIterator(int stream_id)
{
	lock_guard<mutex> lck(mtx);

	m_streams[stream_id].cur_id = 0;

	return true;
}

// ******************************************************************************
size_t CArchive::signature(vector<uint8_t>& v_data)
{
	size_t h = 0;

	for (size_t i = 0; i < v_data.size(); i += 8)
	{
		size_t x = 0;
		for (size_t j = 0; j < 8u && i + j < v_data.size(); ++j)
			x = (x << 8) + (size_t)v_data[i + j];

		x ^= x >> 33;
		x *= 0xff51afd7ed558ccdL;
		x ^= x >> 33;
		x *= 0xc4ceb9fe1a85ec53L;
		x ^= x >> 33;

		h ^= x;
	}

	return h;
}

// ******************************************************************************
bool CArchive::LinkStream(int stream_id, string stream_name, int target_id)
{
	m_streams[stream_id] = m_streams[target_id];
	m_streams[stream_id].stream_name = stream_name;

	return true;
}

// EOF
