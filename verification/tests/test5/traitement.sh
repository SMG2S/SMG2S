#!/bin/bash

grep "@> The eigenvalue" vector1_results.txt > tmp.txt

awk '{print $5 " " $7 }' tmp.txt > tmp2.txt

awk '{print substr($0, 1, length($0)-1)}' tmp2.txt > vector1_results_clean.txt

rm tmp.txt tmp2.txt

tail -n+3 vector1.txt > tmp2.txt

awk '{print $2 " " $3}' tmp2.txt > vector1_initial_clean.txt

rm tmp2.txt



