nThreads=$(grep processor /proc/cpuinfo | wc -l); # Get number of threads
let nParallel=$nThreads*3; # we want to use 0.5 of all available cores max, assuming nThreads is nThreads (not nCores)
let nParallel=$nParallel+7; # The + 7 let's truncating arithmetic round up
let nParallel=$nParallel/8;

MatlabPath=$(which matlab); # find Matlab path
DataParPath=/mnt/s4e_data/RAD/share/EPAD500_new/analysis/DataPar.json; # CHANGE DIR
cd $xASLdir
cd CustomScripts/EPAD

for (( i=1; i<=$nParallel; i++ ));
#do echo Running instance $i";
#here we only run the first (==structural) module (last 1 input argument)
do screen -dmSL EPAD_QC$i nice -n 10 $MatlabPath -nodesktop -nosplash -r "addpath(pwd);cd ../..;ExploreASL_Master_EPAD('$DataParPath', true, true, $i, $nParallel)" &
done
