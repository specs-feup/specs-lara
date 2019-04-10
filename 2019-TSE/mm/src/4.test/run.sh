#!/usr/bin/env bash

runs=1

rm time.csv energy.csv

for m_size in 128 256 512 1024 2018
do
	for t_size in 0 16 32 64 128 256 512
	do
		./mm ${runs} ${m_size} ${t_size} 2>> energy.csv 1>> time.csv
	done
	echo "" >> energy.csv
	echo "" >> time.csv
done
