#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

find . -regex '.*\.PLASMA_VS_TUMOR_RESULT.csv' | grep -v _dis | xargs md5sum 
find . -regex '.*\.PLASMA_VS_TUMOR_RESULT.csv.HBCs.txt' | grep -v _dis | xargs md5sum 

