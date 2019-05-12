#!/bin/bash

echo "===complex double test==="

yhrun -n 48 -N 2 ./smg2s.exe -SIZE 1000000 -L 10 -C 7 > s48.txt
yhrun -n 96 -N 4 ./smg2s.exe -SIZE 2000000 -L 10 -C 7 > s96.txt
yhrun -n 192 -N 8 ./smg2s.exe -SIZE 4000000 -L 10 -C 7 > s192.txt
yhrun -n 384 -N 16 ./smg2s.exe -SIZE 8000000 -L 10 -C 7 > s384.txt
yhrun -n 768 -N 32 ./smg2s.exe -SIZE 16000000 -L 10 -C 7 > s768.txt
#yhrun -n 1536 -N 64 ./smg2s.exe -SIZE 32000000 -L 10 -C 7 > s1536.txt

yhrun -n 48 -N 2 ./smg2s.exe -SIZE 16000000 -L 10 -C 7 > w48.txt
yhrun -n 96 -N 4 ./smg2s.exe -SIZE 16000000 -L 10 -C 7 > w96.txt
yhrun -n 192 -N 8 ./smg2s.exe -SIZE 16000000 -L 10 -C 7 > w192.txt
yhrun -n 384 -N 16 ./smg2s.exe -SIZE 16000000 -L 10 -C 7 > w384.txt
yhrun -n 768 -N 32 ./smg2s.exe -SIZE 16000000 -L 10 -C 7 > w768.txt
#yhrun -n 1536 -N 64 ./smg2s.exe -SIZE 16000000 -L 10 -C 7 > w1536.txt

echo "===Done==="

