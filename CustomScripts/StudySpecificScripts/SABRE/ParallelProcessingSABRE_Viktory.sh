nThreads=$(grep processor /proc/cpuinfo | wc -l); # Get number of threads
let nParallel=$nThreads*3; # we want to use 0.5 of all available cores max, assuming nThreads is nThreads (not nCores)
let nParallel=$nParallel+7; # The + 7 let's truncating arithmetic round up
let nParallel=$nParallel/8;
let nParallel=3;

MatlabPath=/usr/local/apps/matlab/R2018b/bin/matlab; # find Matlab path
DataParPath=/radshare/SABRE/analysis/DataParameters_HiQ.json; # CHANGE DIR
xASLdir=/radshare/EPAD/scripts/ExploreASL;
cd $xASLdir
let iModule=1; #run first only structural module

for (( i=1; i<=$nParallel; i++ ));


screen -dmSL SABRE$i $MatlabPath -nodesktop -nosplash -r "cd('$xASLdir');ExploreASL_Master('$DataParPath', true, true, $i, $nParallel, $iModule);" &


done