#!/bin/bash

echo "Launching Julia files"
nameFile="COMPARE_SOLVERS_YES_GUESS"
nCores=32
nThreads=32
dateString=$(date +'%d-%m-%Y_%H.%M.%S')

nohup /usr/bin/time -v julia -p $nCores -t $nThreads "Metalens OPTIM - Launcher.jl" > "Outputs/out_$nameFile""_n$nCores""_t$nThreads""_$dateString.out" &
echo "Launching completed"