#!/bin/bash
rm -rf timing_scripts/local_tmp
mkdir timing_scripts/local_tmp

. timing_scripts/set_iters.sh

for ((i=0; i < ITERS; i++))
do
./build/release/Commitsnd -e 1 -n 500 -ip 10.11.100.216 -p 28001
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Commitsnd -e 2 -n 1000 -ip 10.11.100.216 -p 28001
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Commitsnd -e 4 -n 15000 -ip 10.11.100.216 -p 28001
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Commitsnd -e 8 -n 50000 -ip 10.11.100.216 -p 28001
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Commitsnd -e 20 -n 500000 -ip 10.11.100.216 -p 28001
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Commitsnd -e 200 -n 10000000 -ip 10.11.100.216 -p 28001
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Commitsnd -e 400 -n 200000000 -ip 10.11.100.216 -p 28001
sleep 2
done