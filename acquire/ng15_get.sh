#!/bin/bash
UA="beacon-search/1.0 (rtg@dxtra.com; academic pulsar timing study)"
cd "${BEACON_DATA:-$HOME}"
U="https://zenodo.org/api/records/16051178/files/NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz/content"
echo "downloading $(date -u +%T)"
curl -s -L -A "$UA" -o ng15.tar.gz "$U"
echo "got $(stat -c%s ng15.tar.gz) bytes"
file ng15.tar.gz | cut -c1-60
mkdir -p ng15_release
echo "extracting the noise files only"
tar -xzf ng15.tar.gz -C ng15_release --wildcards "*/noise/*" 2>/dev/null || \
  tar -xzf ng15.tar.gz -C ng15_release
echo "noise files: $(find ng15_release -path "*noise*" -name "*.txt" | wc -l)"
find ng15_release -path "*noise*" -name "*pars.txt" | head -3
echo DONE
