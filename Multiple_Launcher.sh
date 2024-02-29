#!/bin/bash

echo "Launching Julia files"
nameFile="TEST"
nCores=16
nThreads=$nCores
dateString=$(date +'%d-%m-%Y_%H.%M')

nohup /usr/bin/time -v julia -p $nCores -t $nThreads "Metalens OPTIM - Launcher.jl" > "Outputs/out_$nameFile""_n$nCores""_$dateString.out" &
echo "Launching completed"