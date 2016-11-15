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
./build/release/Tinyconst -n 1 -c sha-256 -e 32,1,1 -ip $1
sleep 2
done