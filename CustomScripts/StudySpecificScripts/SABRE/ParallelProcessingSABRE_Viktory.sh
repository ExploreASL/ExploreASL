let nParallel=3;

MatlabPath=/usr/local/apps/matlab/R2018b/bin/matlab; # find Matlab path
DataParPath=/radshare/SABRE/analysis/DataParameters_HiQ.json; # CHANGE DIR
xASLdir=/radshare/EPAD/scripts/ExploreASL;
cd $xASLdir
let iModule=1; #run first only structural module

for (( i=1; i<=$nParallel; i++ ));


screen -dmSL SABRE$i nice -n 10 $MatlabPath -nodesktop -nosplash -r "cd('$xASLdir');ExploreASL_Master('$DataParPath', true, true, $i, $nParallel, $iModule);" &


done