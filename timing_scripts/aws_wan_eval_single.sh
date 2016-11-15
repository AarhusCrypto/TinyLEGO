#!/bin/bash
rm -rf timing_scripts/local_tmp
mkdir timing_scripts/local_tmp
echo ------------------------------------------------------------------------------------------
echo ----------------------------------------Sequential----------------------------------------
echo ------------------------------------------------------------------------------------------
. timing_scripts/set_iters.sh
for ((i=0; i < SINGLE_ITERS; i++))
do
./build/release/Tinyeval -n 1 -c aes -e 2,1,1 -ip $1 -t 1 >> timing_scripts/local_tmp/1_seq_aes.txt
sleep 2
done

for ((i=0; i < SINGLE_ITERS; i++))
do
./build/release/Tinyeval -n 1 -c sha-256 -e 6,1,1 -ip $1 -t 1 >> timing_scripts/local_tmp/1_seq_sha256.txt
sleep 2
done

awk '{setup+=$1; pre+=$2; off+=$3; onl+=$4};END{printf("1_seq_aes %.2f & %.2f & %.2f & %.2f \n", setup/NR, pre/NR, off/NR, onl/NR)}' timing_scripts/local_tmp/1_seq_aes.txt
awk '{setup+=$1; pre+=$2; off+=$3; onl+=$4};END{printf("1_seq_sha256 %.2f & %.2f & %.2f & %.2f \n", setup/NR, pre/NR, off/NR, onl/NR)}' timing_scripts/local_tmp/1_seq_sha256.txt