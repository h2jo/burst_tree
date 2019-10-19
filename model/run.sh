#!/bin/bash
gcc main.c -o main.x -lm -O3

folder="./"
n=500000
alpha=1.80
a=3.0
h=100.0
s=4.0

./main.x $folder $n $alpha $a $h $s
