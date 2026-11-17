#!/bin/bash
# 11 nov 2026: descarga de las 10 bibliotecas del plan congelado (commit e353b80), con MD5 de ENA.
cd /data/validation_independent_mocks_20261111
echo -e "set\trun\tread_count\tfile\tmd5_ena\tmd5_local\tstatus" > manifest_download.tsv
while read set acc; do
  info=$(curl -s -4 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$acc&result=read_run&fields=read_count,fastq_ftp,fastq_md5&format=tsv" | tail -1)
  rc=$(echo "$info" | cut -f2); urls=$(echo "$info" | cut -f3); md5s=$(echo "$info" | cut -f4)
  for i in 1 2; do
    u=$(echo "$urls" | cut -d";" -f$i); m=$(echo "$md5s" | cut -d";" -f$i); f=raw/$(basename "$u")
    for try in 1 2 3; do
      curl -s -4 --retry 5 -o "$f" "https://$u"; l=$(md5sum "$f" | cut -d" " -f1)
      [ "$l" = "$m" ] && break
    done
    st=OK; [ "$l" = "$m" ] || st=MD5_FAIL
    echo -e "$set\t$acc\t$rc\t$(basename $f)\t$m\t$l\t$st" >> manifest_download.tsv
  done
done <<LIST
A SRR27790137
A SRR27790136
D ERR14820732
D ERR14820828
D ERR14820924
D ERR14821020
D ERR14821116
E SRR16997123
E SRR16997307
E SRR16997311
LIST
echo "DONE $(date -u)" >> manifest_download.tsv
