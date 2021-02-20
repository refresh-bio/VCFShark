// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>
#include <list>
#include <unordered_map>
#include <nmmintrin.h>
#include <chrono>

#include "params.h"
#include "application.h"

#include "rc.h"
#include "sub_rc.h"
#include "io.h"
#include "utils.h"

using namespace std;
using namespace std::chrono;

CParams params;
CApplication *app;

bool parse_params(int argc, char **argv);
void usage_main();
void usage_compress();
void usage_decompress();

// ******************************************************************************
void usage_main()
{
	cerr << "VCFShark v. 1.1 (2021-02-18)\n";
	cerr << "Usage:\n";
	cerr << "  vcfshark <mode>\n";
	cerr << "Parameters:\n";
	cerr << "  mode - one of:\n";
	cerr << "    compress   - compress VCF file\n";
	cerr << "    decompress - decompress VCF file\n";
}

// ******************************************************************************
void usage_compress()
{
	cerr << "VCFShark v. 1.1 (2021-02-18)\n";
	cerr << "Usage:\n";
	cerr << "  vcfshark compress [options] <input_vcf> <archive>\n";
	cerr << "Parameters:\n";
	cerr << "  input_vcf - path to input VCF (or VCF.GZ or BCF) file\n";
	cerr << "  archive - path to output compressed VCF file\n";
	cerr << "Options:\n";
    cerr << "  -nl <value> - ignore rare variants; value is a limit of alternative alleles (default: " << params.neglect_limit << ")\n";
    cerr << "  -t <value>  - max. no. of compressing threads (default: " << params.no_threads << ")\n";
    cerr << "  -c <value>  - compression level [1, 2, 3] (default: " << params.vcs_compression_level << ")\n";
}

// ******************************************************************************
void usage_decompress()
{
	cerr << "VCFShark v. 1.1 (2021-02-18)\n";
	cerr << "Usage:\n";
	cerr << "  vcfshark decompress [options] <archive> <output_vcf>\n";
	cerr << "Parameters:\n";
	cerr << "  archive   - path to input file with compressed VCF file\n";
	cerr << "  output_vcf - path to output VCF file\n";
    cerr << "Options:\n";
    cerr << "  -b - output BCF file (VCF file by default)\n";
    cerr << "  -c [0-9]   set level of compression of the output bcf (number from 0 to 9; 1 by default; 0 means no compression)\n";
	cerr << "  -t <value>  - max. no. of compressing threads (default: " << params.no_threads << ")\n";
}

// ******************************************************************************
bool parse_params(int argc, char **argv)
{
	if (argc < 2)
	{
		usage_main();

		return false;
	}

	if (string(argv[1]) == "compress")
		params.work_mode = work_mode_t::compress;
	else if (string(argv[1]) == "decompress")
		params.work_mode = work_mode_t::decompress;

	// Compress
	if (params.work_mode == work_mode_t::compress)
	{
		if (argc < 4)
		{
			usage_compress();
			return false;
		}

		int i = 2;
		while (i < argc - 2)
		{
			if (string(argv[i]) == "-nl" && i + 1 < argc - 2)
			{
				params.neglect_limit = atoi(argv[i + 1]);
				i += 2;
			}
			else if (string(argv[i]) == "-t" && i + 1 < argc - 2)
			{
				params.no_threads = atoi(argv[i + 1]);
				i += 2;
			}
			else if (string(argv[i]) == "-c" && i + 1 < argc - 2)
			{
				params.vcs_compression_level = atoi(argv[i + 1]);
				if (params.vcs_compression_level < 1 || params.vcs_compression_level > 3)
					params.vcs_compression_level = 3;
				i += 2;
			}
        }

		params.vcf_file_name = string(argv[i]);
		params.db_file_name = string(argv[i+1]);
	}
	else if (params.work_mode == work_mode_t::decompress)
	{
		if (argc < 4)
		{
			usage_decompress();
			return false;
		}

        int i = 2;
        while (i < argc - 2)
        {
            if (string(argv[i]) == "-b")
            {
                params.out_type = file_type::BCF;
                i++;
            }
			else if (string(argv[i]) == "-t" && i + 1 < argc - 2)
			{
				params.no_threads = atoi(argv[i + 1]);
				i += 2;
			}
			else if (string(argv[i]) == "-c")
            {
                i++;
                if(i >= argc - 2)
                {
                    usage_decompress();
                    return false;
                }
                int tmp = atoi(argv[i]);
                if(tmp < 0 || tmp > 9)
                {
                    usage_decompress();
                    return false;
                }
                else
                {
                    if(tmp)
                        params.bcf_compression_level = argv[i][0];
                    else
                        params.bcf_compression_level = 'u';
                }
                i++;
            }
            else
            {
                cerr << "Unknown option : " << argv[i] << endl;
                usage_decompress();
                return false;
            }         
        }

		params.db_file_name = string(argv[i]);
		params.vcf_file_name = string(argv[i+1]);
	}
	else
	{
		cerr << "Unknown mode : " << argv[2] << endl;
		usage_main();

		return false;
	}

	return true;
}

// ******************************************************************************
int main(int argc, char **argv)
{
	if (!parse_params(argc, argv))
		return 0;

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	app = new CApplication(params);

	bool result = true;

	if (params.work_mode == work_mode_t::compress)
		result = app->CompressDB();
	else if (params.work_mode == work_mode_t::decompress)
		result = app->DecompressDB();

	delete app;

	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

	if (!result)
		std::cout << "Critical error!\n";

	std::cout << "Processing time: " << time_span.count() << " seconds.\n";

	fflush(stdout);

	return 0;
}

// EOF
