# VCFShark:  how to squeeze a VCF file

VCFShark is a tool to compress any VCF file. It achieves compression ratios up to an order of magnitude better than the de facto standards (gzipped VCF and BCF).

As an input it takes a VCF (or VCF.GZ or BCF) file. 

Requirements
--------------

VCFShark requires:

* A modern, C++14 ready compiler such as `g++` version 7.2 or higher or `clang` version 3.4 or higher.
* A 64-bit operating system. Either Mac OS X or Linux are currently supported.

Installation
--------------

To download, build and install VCFShark use the following commands.
```sh
git clone https://github.com/refresh-bio/VCFShark.git
cd VCFShark
./install.sh 
make
make install
```
The `install.sh` script downloads and installs the HTSlib library into the `include` and `lib` directories in the `VCFShark/htslib` directory. The compiled libbsc library is statically linked (.

By default VCFShark is installed in the `bin` directory of the `/usr/local/` directory. A different location prefix can be specified with `prefix` parameter:
```sh
make prefix=/usr/local install
```
---
To uninstall VCFShark:
```sh
make uninstall
```
This uninstalls VCFShark from the `/usr/local` directory. To uninstall from different location use the `prefix` parameter:
```sh
make prefix=/usr/local uninstall
```
To uninstall the HTSlib library use the provided uninstall script:
```sh
./uninstall.sh 
```
---
To clean the VCFShark build use:
```sh
make clean
```

Usage
--------------
* Compress the input VCF/BCF/VCF.gz file.
```
Input: <input_vcf> VCF/BCF file. 
Output: <archive> output file with the compressed VCF

Usage: 
vcfshark compress [options] <input_vcf> <output_db>
Parameters:
  input_vcf - path to input VCF (or VCF.GZ or BCF) file
  archive - path to the output compressed VCF
Options:
  -nl <value> - ignore rare variants; value is a limit of alternative alleles (default: 10)
  -t <value>  - max. no. of compressing threads (default: 8)
  ```
  
 * Decompress the archive.
 ```
Input: <archive> archive 
Output: <output_vcf> VCF/BCF file.
 
Usage: 
vcfshark decompress [options] <archive> <output_vcf>
Parameters:
  archive   - path to compressed VCF
  output_vcf - path to output VCF/BCF file
Options:
  -b - output BCF file (VCF file by default)
  -c [0-9]   set level of compression of the output bcf (number from 0 to 9; 1 by default; 0 means no compression)	
  -t <value>  - max. no. of compressing threads (default: 8)
 ```
 
 
Toy example
--------------

There is an example VCF files in the `toy_ex` folder: `toy.vcf`. It can be used to test VCFShark. 
The instructions should be called within `toy_ex` folder.

To compress a single example VCF file and store the archive called `toy.vcfshark` in the `toy_ex` folder:
```sh
../vcfshark compress toy.vcf toy.vcfshark
```
This will create an archive `toy.vcfshark`.

To decompress the `toy.vcfshark`  archive to a VCF file `toy_decomp.vcf`:
```sh
../vcfshark decompress toy.vcfshark toy_decomp.vcf
```

For more options see Usage section.


Dockerfile
--------------
Dockerfile can be used to build a Docker image with all necessary dependencies and VCFShark compressor. 

The first image is based on Ubuntu 18.04 (Dockerfile_ubuntu), the second one on CentOS 7 (Dockerfile_centos). 

To build a Docker image and run a Docker container, you need Docker Desktop (https://www.docker.com). 

Example commands (run it within a directory with Dockerfile):
```sh
docker build -f Dockerfile_ubuntu -t ubuntu-vcfshark .
docker run -it ubuntu-vcfshark
```
or:
```sh
docker build -f Dockerfile_centos -t centos-vcfshark .
docker run -it centos-vcfshark
```

Note: The Docker image is not intended as a way of using VCFShark. It can be used to test the instalation process and basic operation of the VCFShark tool.



Developers
--------------
The VCFShark algorithm was invented by [Sebastian Deorowicz](https://github.com/sebastiandeorowicz) and [Agnieszka Danek](https://github.com/agnieszkadanek).
The implementation is by Sebastian Deorowicz and Agnieszka Danek.



