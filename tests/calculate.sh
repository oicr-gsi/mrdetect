#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

find -name SNP.count.txt | grep -v _dis | xargs md5sum 
find -name *.HBCs.csv | grep -v _dis | xargs md5sum 
find -name *.mrdetect.vaf.txt | grep -v _dis | xargs md5sum
find -name *.mrdetect.txt | grep -v _dis | xargs md5sum
