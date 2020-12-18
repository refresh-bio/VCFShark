#!/bin/sh

INSTALL_DIR=`pwd`/htslib


#if [ ! -z "${1}" ]; then
#   INSTALL_DIR=${1}
#fi

     
case "$(uname -s)" in

   Darwin)

     mkdir htslib
     cd htslib
     wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
     tar -xf htslib-1.10.2.tar.bz2
     cd htslib-1.10.2
     ./configure --libdir=${INSTALL_DIR}/lib CC=clang
     make 
     make prefix=${INSTALL_DIR} install
     cd ..
     cd ..
     ;;

   Linux)
   
     mkdir htslib
     cd htslib
     wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
     tar -xf htslib-1.10.2.tar.bz2
     cd htslib-1.10.2
     ./configure --libdir=${INSTALL_DIR}/lib  
     make
     make prefix=${INSTALL_DIR} install
     cd ..
     cd ..
     ;;


   *)
     echo 'unsupported OS' 
     ;;
esac
