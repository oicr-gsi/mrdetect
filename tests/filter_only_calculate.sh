#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

find -name SNP.count.txt | grep -v _dis | xargs md5sum 
find . -name *.vcf | grep -v ^# | xargs md5sum
