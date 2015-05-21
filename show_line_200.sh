#!/bin/bash

for f in all_snp*april.txt
do
	echo "$f"
	sed '200q;d' $f
done
