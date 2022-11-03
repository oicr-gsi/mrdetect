#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

find . -regex '.*\.json.*' | grep -v _dis | xargs md5sum 

