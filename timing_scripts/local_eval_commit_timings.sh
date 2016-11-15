#!/bin/bash
rm -rf timing_scripts/local_tmp
mkdir timing_scripts/local_tmp

. timing_scripts/set_iters.sh

for ((i=0; i < ITERS; i++))
do
./build/release/Commitrec -e 1 -n 500 -ip 10.11.100.216 -p 28001 -t 1 >> timing_scripts/local_tmp/500_commit.txt
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Commitrec -e 2 -n 1000 -ip 10.11.100.216 -p 28001 -t 1 >> timing_scripts/local_tmp/1000_commit.txt
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Commitrec -e 4 -n 15000 -ip 10.11.100.216 -p 28001 -t 1 >> timing_scripts/local_tmp/15000_commit.txt
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Commitrec -e 8 -n 50000 -ip 10.11.100.216 -p 28001 -t 1 >> timing_scripts/local_tmp/50000_commit.txt
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Commitrec -e 20 -n 500000 -ip 10.11.100.216 -p 28001 -t 1 >> timing_scripts/local_tmp/500000_commit.txt
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Commitrec -e 200 -n 10000000 -ip 10.11.100.216 -p 28001 -t 1 >> timing_scripts/local_tmp/10000000_commit.txt
sleep 2
done

for ((i=0; i < ITERS; i++))
do
./build/release/Commitrec -e 400 -n 200000000 -ip 10.11.100.216 -p 28001 -t 1 >> timing_scripts/local_tmp/200000000_commit.txt
sleep 2
done

awk '{num_commit+=$1; num_thread+=$2; commit+=$3; commit_ot+=$4; decommit+=$5};END{printf("500_commit %d & %d & %.2f & (%.2f) & %.2f \n", num_commit/NR, num_thread/NR, commit/NR, commit_ot/NR, decommit/NR)}' timing_scripts/local_tmp/500_commit.txt

awk '{num_commit+=$1; num_thread+=$2; commit+=$3; commit_ot+=$4; decommit+=$5};END{printf("1000_commit %d & %d & %.2f & (%.2f) & %.2f \n", num_commit/NR, num_thread/NR, commit/NR, commit_ot/NR, decommit/NR)}' timing_scripts/local_tmp/1000_commit.txt

awk '{num_commit+=$1; num_thread+=$2; commit+=$3; commit_ot+=$4; decommit+=$5};END{printf("15000_commit %d & %d & %.2f & (%.2f) & %.2f \n", num_commit/NR, num_thread/NR, commit/NR, commit_ot/NR, decommit/NR)}' timing_scripts/local_tmp/15000_commit.txt

awk '{num_commit+=$1; num_thread+=$2; commit+=$3; commit_ot+=$4; decommit+=$5};END{printf("50000_commit %d & %d & %.2f & (%.2f) & %.2f \n", num_commit/NR, num_thread/NR, commit/NR, commit_ot/NR, decommit/NR)}' timing_scripts/local_tmp/50000_commit.txt

awk '{num_commit+=$1; num_thread+=$2; commit+=$3; commit_ot+=$4; decommit+=$5};END{printf("500000_commit %d & %d & %.2f & (%.2f) & %.2f \n", num_commit/NR, num_thread/NR, commit/NR, commit_ot/NR, decommit/NR)}' timing_scripts/local_tmp/500000_commit.txt

awk '{num_commit+=$1; num_thread+=$2; commit+=$3; commit_ot+=$4; decommit+=$5};END{printf("10000000_commit %d & %d & %.2f & (%.2f) & %.2f \n", num_commit/NR, num_thread/NR, commit/NR, commit_ot/NR, decommit/NR)}' timing_scripts/local_tmp/10000000_commit.txt

awk '{num_commit+=$1; num_thread+=$2; commit+=$3; commit_ot+=$4; decommit+=$5};END{printf("200000000_commit %d & %d & %.2f & (%.2f) & %.2f \n", num_commit/NR, num_thread/NR, commit/NR, commit_ot/NR, decommit/NR)}' timing_scripts/local_tmp/200000000_commit.txt