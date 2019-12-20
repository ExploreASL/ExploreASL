Threads=$(grep processor /proc/cpuinfo | wc -l); # Get number of threads
let nParallel=16;

MatlabPath=$(which matlab); # find Matlab path
DataParPath=/home/hjmutsaerts/scratch/Data/DataPar.json; # CHANGE DIR
cd $xASLdir
cd CustomScripts/EPAD

for (( i=1; i<=$nParallel; i++ ));
#do echo Running instance $i";
#here we only run the first (==structural) module (last 1 input argument)
do screen -dmSL EPAD_QC$i nice -n 10 $MatlabPath -nodesktop -nosplash -r "addpath(pwd);cd ../..;ExploreASL_Master_EPAD('$DataParPath', true, true, $i, $nParallel, [3])" &
done
