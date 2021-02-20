## The scripts in this folder can be used to run VCFShark on real data. 

### Each folder contains a script to:
1. Download a VCF file, 
2. If needed: repair the downloaded VCF file (so it complies with the VCF specification) using bcftools.
3. Compress the VCF file using VCFShark (input format: vcf.gz).
4. Decompress the VCF file using VCFShark (output format: vcf.gz).
5. Running "check" script on input and output files, to check if the compression/decompression works. 

--------------

### Example procedure to run VCFShark on 1000GPp1 chr22 VCF file and on Tomato VCF file. 

(Note: on macOS it may be necessary to change `MD5` variable in `check` script to `md5`.)

1. Build VCFShark in the main VCFShark folder, according to instructions:
```
./install.sh 
make
```
2. Change directory to `large-ex`, download and build bcftools (required by `check` script and to fix some of the VCF files).
```
cd large-ex
./get_bcftools
```
3. Change the directory to folder with chosen data and run the available script:
```
cd 1000GPp1
./run_1000GPp1
```
4. Next, to run VCFShark on Tomato dataset:
```
cd ../Tomato
./run_Tomato
```

--------------

### Comments:
- The script `get_bcftools` downloads and install the bcftools in the current directory (`large_ex` folder). 
- It is necessary to download and build bcftools, and set BCFTOOLS variable (in `run_*` and `check` scripts) accordingly (by defaut is is assumed it is in the `large_ex` folder). 
- It is necessary to set VCFSHARK variable (in `run_*` scripts). By default the VCFShark variable in all scripts refers to the executable in the main `VCFShark` folder.
- As for VCF repairing: some VCF files does not contain the meta-information lines in the header, which are recommended in the VCF specification and required by the binary counterpart of VCF (BCF) and by VCFShark. (From VCF specification: "Note that BCF, the binary counterpart of VCF, requires that all entries are present.")
Bcftools is used to fix this issue, adding missing header lines to the input VCF file. Missing header lines, if required, are in `missing_header` file in related folder.
- As for checking: the order of subfields in columns of the VCF file compressed/decompressed with VCFhark match the order of the occurrence of meta-information lines related to that subfields in the VCF header.
Therefore the order of subfields may differ in the input and output VCF files. Note that it has no influence on the data contained in the VCF file. 
Also, VCFShark, when outputting the VCF file from the archive, cuts the isignificant zeros (for exmple: `8.63300e+06` becomes `8.633e+06`).  
   Taking into account the described issues, to compare the input and output file:
the script "check" uses bcftools and gawk to create two VCF files (text files): in.vcf and out.vcf, based on input and output VCF/BCF/VCF.GZ files, respectively. 
The only modifications are: sorting the subfields in the INFO column and removing insignificant zeros from the input file.
Next, the md5 of both VCF files are calculated and compared. 
- The utility to check md5 of a file in the `check` script is set to `md5sum` by default; to change it, it is necessary to change MD5 variable in `check` script.


--------------

### Disk space usage

The scripts do not remove any created files. Running full tests with the scripts provided require large amount of space. \
See the approximate space requirements for full test on each VCF file:

#### 1000GPp1
Total: 
- 28 GB  

Including:
- 1.8 GB   ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz
- 0.6 GB   1000GPp1_chr22.vcfshark
- 11 GB    in.vcf
- 11 GB    out.vcf


#### 1000GPp3
Total: 
- 22 GB  

Including:
- 0.2 GB   ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
- 0.02 GB   1000GPp3_chr22.vcfshark
- 11 GB    in.vcf
- 11 GB    out.vcf

#### ExAC
Total: 
- 97 GB  

Including:
- 4.6 GB   ExAC.r1.sites.vep.vcf.gz 
- 1.9 GB   exac.vcfshark
- 38 GB    in.vcf
- 38 GB    out.vcf

#### Tomato
Total: 
- 26 GB  

Including:
- 1.0 GB  360_merged_2.50.vcf.gz
- 0.4 GB   tomato.vcfshark
- 12 GB    in.vcf
- 12 GB    out.vcf

#### Rice
Total:   
- 56 GB  

Including:
- 2.0 GB   IRIS_313-10000.snp.vcf.gz
- 0.4 GB   rice.vcfshark
- 26 GB    in.vcf
- 26 GB    out.vcf

#### gnomAD2
Total:  
- 108 GB  

Including:
- 6.5 GB   gnomad.genomes.r2.1.1.sites.22.vcf.bgz
- 0.9 GB   gnomad2.vcfshark
- 46 GB    in.vcf
- 46 GB    out.vcf
