nThreads=$(grep processor /proc/cpuinfo | wc -l); # Get number of threads
let nParallel=$nThreads*3; # we want to use 0.5 of all available cores max, assuming nThreads is nThreads (not nCores)
let nParallel=$nParallel+7; # The + 7 let's truncating arithmetic round up
let nParallel=$nParallel/8;

MatlabPath=$(which matlab); # find Matlab path
DataParPath=/Users/henk/surfdrive/SABRE/analysis/DataParameters_HiQ.json; # CHANGE DIR
cd $xASLdir
cd CustomScripts/SABRE

for (( i=1; i<=$nParallel; i++ ));
#do echo Running instance $i";
do screen -dmSL SABRE$i nice -n 10 $MatlabPath -nodesktop -nosplash -r "addpath(pwd);cd ../..;ExploreASL_Master('$DataParPath', true, true, $i, $nParallel);" &
done
