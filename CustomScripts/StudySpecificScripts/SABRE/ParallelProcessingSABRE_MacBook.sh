nThreads=$(grep processor /proc/cpuinfo | wc -l); # Get number of threads
let nParallel=$nThreads*3; # we want to use 0.5 of all available cores max, assuming nThreads is nThreads (not nCores)
let nParallel=$nParallel+7; # The + 7 let's truncating arithmetic round up
let nParallel=$nParallel/8;
let nParallel=8;

MatlabPath=/Applications/MATLAB_R2016b.app/bin/matlab; # find Matlab path
DataParPath=/Users/henk/ExploreASL/ASL/SABRE/Analysis2/DataParameters_HiQ.json; # CHANGE DIR
xASLdir=/Users/henk/ExploreASL/ExploreASL;
cd $xASLdir

for (( i=1; i<=$nParallel; i++ ));
do echo Running instance $i";
#do screen -dmSL SABRE$i $MatlabPath -nodesktop -nosplash -r "ExploreASL_Master('$DataParPath', true, true, $i, $nParallel);" &
done

let iModule=2; # run ASL module, skip DARTEL module

let nParallel=8;

i=1;
screen -dmSL SABRE$i $MatlabPath -nodesktop -nosplash -r "cd $xASLdir;ExploreASL_Master('$DataParPath', true, true, $i, $nParallel, $iModule);" &
i=2;
screen -dmSL SABRE$i $MatlabPath -nodesktop -nosplash -r "cd $xASLdir;ExploreASL_Master('$DataParPath', true, true, $i, $nParallel, $iModule);" &
i=3;
screen -dmSL SABRE$i $MatlabPath -nodesktop -nosplash -r "cd $xASLdir;ExploreASL_Master('$DataParPath', true, true, $i, $nParallel, $iModule);" &
i=4;
screen -dmSL SABRE$i $MatlabPath -nodesktop -nosplash -r "cd $xASLdir;ExploreASL_Master('$DataParPath', true, true, $i, $nParallel, $iModule);" &
i=5;
screen -dmSL SABRE$i $MatlabPath -nodesktop -nosplash -r "cd $xASLdir;ExploreASL_Master('$DataParPath', true, true, $i, $nParallel, $iModule);" &
i=6;
screen -dmSL SABRE$i $MatlabPath -nodesktop -nosplash -r "cd $xASLdir;ExploreASL_Master('$DataParPath', true, true, $i, $nParallel, $iModule);" &
i=7;
screen -dmSL SABRE$i $MatlabPath -nodesktop -nosplash -r "cd $xASLdir;ExploreASL_Master('$DataParPath', true, true, $i, $nParallel, $iModule);" &
i=8;
screen -dmSL SABRE$i $MatlabPath -nodesktop -nosplash -r "cd $xASLdir;ExploreASL_Master('$DataParPath', true, true, $i, $nParallel, $iModule);" &