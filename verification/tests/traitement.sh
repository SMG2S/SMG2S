#!/bin/bash

grep "@> The eigenvalue" vector1_results.txt > tmp.txt

awk '{print $5 " " $7 }' tmp.txt > tmp2.txt

awk '{print substr($0, 1, length($0)-1)}' tmp2.txt > tmp3.txt

awk '{print NR  " " $0}' tmp3.txt > tmp4.txt

NB=`wc -l tmp4.txt | awk '{print $1}'`

awk 'BEGIN{print '$NB' " " '$NB' " " '$NB'}{print}' tmp4.txt  > tmp5.txt

awk 'BEGIN{print "%%MatrixMarket matrix coordinate real general"}{print}' tmp5.txt   > vector1_results_clean.txt

rm tmp.txt tmp2.txt tmp3.txt tmp4.txt tmp5.txt
