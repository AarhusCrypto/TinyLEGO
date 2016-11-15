#!/bin/bash
echo ------------------------------------------------------------------------------------------
echo ----------------------------------------Sequential----------------------------------------
echo ------------------------------------------------------------------------------------------
. timing_scripts/set_iters.sh
for ((i=0; i < SINGLE_ITERS; i++))
do
./build/release/Tinyconst -n 1 -c aes -e 8,1,1 -ip $1
sleep 2
done

for ((i=0; i < SINGLE_ITERS; i++))
do
./build/release/Tinyconst -n 1 -c sha-256 -e 16,1,1 -ip $1
sleep 2
done

echo ------------------------------------------------------------------------------------------
for ((i=0; i < ITERS; i++))
do
./build/release/Tinyconst -n 32 -c aes -e 72,1,1 -ip $1
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Tinyconst -n 32 -c sha-256 -e 256,1,1 -ip $1
sleep 2
done

echo ------------------------------------------------------------------------------------------
for ((i=0; i < ITERS; i++))
do
./build/release/Tinyconst -n 128 -c aes -e 128,1,1 -ip $1
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Tinyconst -n 128 -c sha-256 -e 1024, 1, 1 -ip $1
sleep 2
done

echo ------------------------------------------------------------------------------------------
for ((i=0; i < ITERS; i++))
do
./build/release/Tinyconst -n 1024 -c aes -e 400,1,1 -ip $1
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Tinyconst -n 256 -c sha-256 -e 1024,1,1 -ip $1
sleep 2
done

echo ------------------------------------------------------------------------------------------
echo -----------------------------------------Parallel-----------------------------------------
echo ------------------------------------------------------------------------------------------

#Parallel
for ((i=0; i < ITERS; i++))
do
./build/release/Tinyconst -n 32 -c aes -e 72,32,32 -ip $1
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Tinyconst -n 32 -c sha-256 -e 256,32,32 -ip $1
sleep 2
done

echo ------------------------------------------------------------------------------------------
for ((i=0; i < ITERS; i++))
do
./build/release/Tinyconst -n 128 -c aes -e 128,128,128 -ip $1
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Tinyconst -n 128 -c sha-256 -e 400,128,128 -ip $1
sleep 2
done

echo ------------------------------------------------------------------------------------------
for ((i=0; i < ITERS; i++))
do
./build/release/Tinyconst -n 1024 -c aes -e 400,400,512 -ip $1
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Tinyconst -n 256 -c sha-256 -e 1024,256,256 -ip $1
sleep 2
done