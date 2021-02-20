#pragma once
// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include <cstdio>
#include <vector>
#include <map>
#include <list>
#include <string>
#include <thread>
#include <mutex>
#include <unordered_map>

using namespace std;

// ************************************************************************************
class CArchive
{
	bool input_mode;
	FILE* f;
	size_t f_offset;
	string file_name;

	struct part_t{
		size_t offset;
		size_t size;

		part_t() : offset(0), size(0)
		{};

		part_t(size_t _offset, size_t _size) : offset(_offset), size(_size)
		{};
	};

	typedef struct {
		string stream_name;
		size_t cur_id;
		size_t raw_size;
		vector<part_t> parts;
		vector<uint64_t> signatures;
	} stream_t;

	map<int, stream_t> m_streams;
	mutex mtx;

	unordered_map<size_t, pair<int, int>> uo_signatures;

	bool serialize();
	bool deserialize();
	size_t write_fixed(size_t x, FILE* file);
	size_t write(size_t x, FILE *file);
	size_t write(string s, FILE* file);
	size_t read_fixed(size_t& x, FILE* file);
	size_t read(size_t& x, FILE* file);
	size_t read(string& s, FILE* file);
	size_t signature(vector<uint8_t>& v_data);

public:
	CArchive(bool _input_mode);
	~CArchive();

	bool Open(string _file_name);
	bool Close();

	int RegisterStream(string stream_name);
	int GetStreamId(string stream_name);

	bool AddPart(int stream_id, vector<uint8_t> &v_data, size_t metadata = 0);
	int AddPartPrepare(int stream_id);
	bool AddPartComplete(int stream_id, int part_id, vector<uint8_t>& v_data, size_t metadata = 0);

	bool GetPart(int stream_id, vector<uint8_t> &v_data, size_t &metadata);
	void SetRawSize(int stream_id, size_t raw_size);
	size_t GetRawSize(int stream_id);
	size_t GetCompressedSize(int stream_id);
	bool ResetStreamPartIterator(int stream_id);

	bool LinkStream(int stream_id, string stream_name, int target_id);


	size_t GetNoStreams()
	{
		lock_guard<mutex> lck(mtx);

		return m_streams.size();
	}
};

// EOF
