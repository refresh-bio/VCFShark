all: vcfshark

VCFShark_ROOT_DIR=.
VCFShark_MAIN_DIR=src
LIBS_DIR=/usr/local/lib
INCLUDE_DIR=libbsc
HTS_INCLUDE_DIR=htslib/include
HTS_LIB_DIR=htslib/lib

CC 	= g++
CFLAGS	= -Wall -O3 -m64 -std=c++14 -pthread -mavx -I $(HTS_INCLUDE_DIR) -I $(INCLUDE_DIR) -fpermissive
CLINK	= -lm -O3 -std=c++14 -pthread -mavx -lz -lbz2 -lcurl -lssl -lcrypto -llzma -L $(LIBS_DIR) 

ifdef MSVC     # Avoid the MingW/Cygwin sections
    uname_S := Windows
else                          # If uname not available => 'not' 
    uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
endif

ifeq ($(uname_S),Linux)
	BSC_LIB_DIR=libbsc/linux
endif

ifeq ($(uname_S),Darwin)
	BSC_LIB_DIR=libbsc/macos
endif

# default install location (binary placed in the /bin folder)
prefix      = /usr/local

# optional install location
exec_prefix = $(prefix)


%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

vcfshark: $(VCFShark_MAIN_DIR)/application.o \
	$(VCFShark_MAIN_DIR)/archive.o \
	$(VCFShark_MAIN_DIR)/bsc.o \
	$(VCFShark_MAIN_DIR)/buffer.o \
	$(VCFShark_MAIN_DIR)/cfile.o \
	$(VCFShark_MAIN_DIR)/cfile_impl.o \
	$(VCFShark_MAIN_DIR)/format.o \
	$(VCFShark_MAIN_DIR)/graph_opt.o \
	$(VCFShark_MAIN_DIR)/main.o \
	$(VCFShark_MAIN_DIR)/pbwt.o \
	$(VCFShark_MAIN_DIR)/text_pp.o \
	$(VCFShark_MAIN_DIR)/utils.o \
	$(VCFShark_MAIN_DIR)/vcf.o
	$(CC) -o $(VCFShark_ROOT_DIR)/$@  \
	$(VCFShark_MAIN_DIR)/application.o \
	$(VCFShark_MAIN_DIR)/archive.o \
	$(VCFShark_MAIN_DIR)/bsc.o \
	$(VCFShark_MAIN_DIR)/buffer.o \
	$(VCFShark_MAIN_DIR)/cfile.o \
	$(VCFShark_MAIN_DIR)/cfile_impl.o \
	$(VCFShark_MAIN_DIR)/format.o \
	$(VCFShark_MAIN_DIR)/graph_opt.o \
	$(VCFShark_MAIN_DIR)/main.o \
	$(VCFShark_MAIN_DIR)/pbwt.o \
	$(VCFShark_MAIN_DIR)/text_pp.o \
	$(VCFShark_MAIN_DIR)/utils.o \
	$(VCFShark_MAIN_DIR)/vcf.o \
	$(BSC_LIB_DIR)/libbsc.a \
	$(HTS_LIB_DIR)/libhts.a \
	$(CLINK)

clean:
	-rm $(VCFShark_MAIN_DIR)/*.o
	-rm vcfshark

install:
	mkdir -p -m 755 $(exec_prefix)/bin
	cp vcfshark $(exec_prefix)/bin/

uninstall:
	rm  $(exec_prefix)/bin/vcfshark


