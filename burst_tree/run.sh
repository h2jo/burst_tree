#!/bin/bash
gcc main.c -o main.x -lm -O3

fn="example"
folder="./"

./main.x $fn $folder
