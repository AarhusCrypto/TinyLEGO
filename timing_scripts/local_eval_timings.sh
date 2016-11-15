#!/bin/bash
rm -rf timing_scripts/local_tmp
mkdir timing_scripts/local_tmp
echo ------------------------------------------------------------------------------------------
echo ----------------------------------------Sequential----------------------------------------
echo ------------------------------------------------------------------------------------------
. timing_scripts/set_iters.sh
for ((i=0; i < SINGLE_ITERS; i++))
do
./build/release/Tinyeval -n 1 -c aes -e 8,1,1 -ip $1 -t 1 >> timing_scripts/local_tmp/1_seq_aes.txt
sleep 2
done

for ((i=0; i < SINGLE_ITERS; i++))
do
./build/release/Tinyeval -n 1 -c sha-256 -e 16,1,1 -ip $1 -t 1 >> timing_scripts/local_tmp/1_seq_sha256.txt
sleep 2
done

echo ------------------------------------------------------------------------------------------
for ((i=0; i < ITERS; i++))
do
./build/release/Tinyeval -n 32 -c aes -e 64,1,1 -ip $1 -t 1 >> timing_scripts/local_tmp/32_seq_aes.txt
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Tinyeval -n 32 -c sha-256 -e 128,1,1 -ip $1 -t 1 >> timing_scripts/local_tmp/32_seq_sha256.txt
sleep 2
done

echo ------------------------------------------------------------------------------------------

for ((i=0; i < ITERS; i++))
do
./build/release/Tinyeval -n 128 -c aes -e 96,1,1 -ip $1 -t 1 >> timing_scripts/local_tmp/128_seq_aes.txt
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Tinyeval -n 128 -c sha-256 -e 128,1,1 -ip $1 -t 1 >> timing_scripts/local_tmp/128_seq_sha256.txt
sleep 2
done

echo ------------------------------------------------------------------------------------------

for ((i=0; i < ITERS; i++))
do
./build/release/Tinyeval -n 1024 -c aes -e 160,1,1 -ip $1 -t 1 >> timing_scripts/local_tmp/1024_seq_aes.txt
sleep 2
done

echo ------------------------------------------------------------------------------------------
echo -----------------------------------------Parallel-----------------------------------------
echo ------------------------------------------------------------------------------------------

for ((i=0; i < ITERS; i++))
do
./build/release/Tinyeval -n 32 -c aes -e 64,32,32 -ip $1 -t 1 >> timing_scripts/local_tmp/32_par_aes.txt
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Tinyeval -n 32 -c sha-256 -e 128,32,32 -ip $1 -t 1 >> timing_scripts/local_tmp/32_par_sha256.txt
sleep 2
done

echo ------------------------------------------------------------------------------------------

for ((i=0; i < ITERS; i++))
do
./build/release/Tinyeval -n 128 -c aes -e 96,128,128 -ip $1 -t 1 >> timing_scripts/local_tmp/128_par_aes.txt
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Tinyeval -n 128 -c sha-256 -e 128,128,128 -ip $1 -t 1 >> timing_scripts/local_tmp/128_par_sha256.txt
sleep 2
done

echo ------------------------------------------------------------------------------------------

for ((i=0; i < ITERS; i++))
do
./build/release/Tinyeval -n 1024 -c aes -e 160,200,128 -ip $1 -t 1 >> timing_scripts/local_tmp/1024_par_aes.txt
sleep 2
done

awk '{setup+=$1; pre+=$2; off+=$3; onl+=$4};END{printf("1_seq_aes %.2f & %.2f & %.2f & %.2f \n", setup/NR, pre/NR, off/NR, onl/NR)}' timing_scripts/local_tmp/1_seq_aes.txt
awk '{setup+=$1; pre+=$2; off+=$3; onl+=$4};END{printf("1_seq_sha256 %.2f & %.2f & %.2f & %.2f \n", setup/NR, pre/NR, off/NR, onl/NR)}' timing_scripts/local_tmp/1_seq_sha256.txt

awk '{setup+=$1; pre+=$2; off+=$3; onl+=$4};END{printf("32_seq_aes %.2f & %.2f & %.2f & %.2f \n", setup/NR, pre/NR, off/NR, onl/NR)}' timing_scripts/local_tmp/32_seq_aes.txt
awk '{setup+=$1; pre+=$2; off+=$3; onl+=$4};END{printf("32_seq_sha256 %.2f & %.2f & %.2f & %.2f \n", setup/NR, pre/NR, off/NR, onl/NR)}' timing_scripts/local_tmp/32_seq_sha256.txt

awk '{setup+=$1; pre+=$2; off+=$3; onl+=$4};END{printf("128_seq_aes %.2f & %.2f & %.2f & %.2f \n", setup/NR, pre/NR, off/NR, onl/NR)}' timing_scripts/local_tmp/128_seq_aes.txt
awk '{setup+=$1; pre+=$2; off+=$3; onl+=$4};END{printf("128_seq_sha256 %.2f & %.2f & %.2f & %.2f \n", setup/NR, pre/NR, off/NR, onl/NR)}' timing_scripts/local_tmp/128_seq_sha256.txt

awk '{setup+=$1; pre+=$2; off+=$3; onl+=$4};END{printf("1024_seq_aes %.2f & %.2f & %.2f & %.2f \n", setup/NR, pre/NR, off/NR, onl/NR)}' timing_scripts/local_tmp/1024_seq_aes.txt

awk '{setup+=$1; pre+=$2; off+=$3; onl+=$4};END{printf("32_par_aes %.2f & %.2f & %.2f & %.2f \n", setup/NR, pre/NR, off/NR, onl/NR)}' timing_scripts/local_tmp/32_par_aes.txt
awk '{setup+=$1; pre+=$2; off+=$3; onl+=$4};END{printf("32_par_sha256 %.2f & %.2f & %.2f & %.2f \n", setup/NR, pre/NR, off/NR, onl/NR)}' timing_scripts/local_tmp/32_par_sha256.txt

awk '{setup+=$1; pre+=$2; off+=$3; onl+=$4};END{printf("128_par_aes %.2f & %.2f & %.2f & %.2f \n", setup/NR, pre/NR, off/NR, onl/NR)}' timing_scripts/local_tmp/128_par_aes.txt
awk '{setup+=$1; pre+=$2; off+=$3; onl+=$4};END{printf("128_par_sha256 %.2f & %.2f & %.2f & %.2f \n", setup/NR, pre/NR, off/NR, onl/NR)}' timing_scripts/local_tmp/128_par_sha256.txt

awk '{setup+=$1; pre+=$2; off+=$3; onl+=$4};END{printf("1024_par_aes %.2f & %.2f & %.2f & %.2f \n", setup/NR, pre/NR, off/NR, onl/NR)}' timing_scripts/local_tmp/1024_par_aes.txt