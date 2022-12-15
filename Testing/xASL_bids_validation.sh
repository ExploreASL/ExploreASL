#!/bin/bash
# Simple script to run a bids-validator on the TestFlavorDatabase 

# input arguments 
checkDir=/home/radv/mshammer/my-scratch/FlavorDatabase/
xASLDir=home/radv/mshammer/ExploreAsl

ARGV=("$@")
ARGC=("$#")

cd ${checkDir}
FolderArray=(*/)
lengthDir="${#FolderArray[@]}"
 
today=$(date +"%FT%H:%M%:z") 
ResultDirToday=${checkDir}/${today}
mkdir ${ResultDirToday}

for ((i=0; i<${lengthDir}; i++)); do
    bids-validator ${checkDir}/${FolderArray[i]}rawdataReference >> ${ResultDirToday}/${FolderArray[i]::-1}.json --no-color --json 
done

#matlab 