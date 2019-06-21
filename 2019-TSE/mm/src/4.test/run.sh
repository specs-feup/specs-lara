#!/usr/bin/env bash

runs=5

rm time.csv energy.csv

for m_size in 512 1024 2048 4096 8192
do
	echo -e "matrix: ${m_size}"
	for t_size in 0 64 128 256 512 1024
	do
		echo -e "   tile: ${t_size}"
		./mm ${runs} ${m_size} ${t_size} 2>> energy.csv 1>> time.csv
	done
	echo "" >> energy.csv
	echo "" >> time.csv
done
