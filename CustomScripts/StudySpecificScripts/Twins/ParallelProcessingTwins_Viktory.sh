let nParallel=10;
MatlabPath=/usr/local/apps/matlab/R2018b/bin/matlab;
DataParPath=/radshare/Twins/Twins_ExploreASL/DATA_PAR_TwinStudy.m;
xASLdir=/radshare/EPAD/scripts/ExploreASL;
cd $xASLdir
for (( i=1; i<=$nParallel; i++ ));
do screen -dmSL TWINS$i nice -n 10 $MatlabPath -nodesktop -nosplash -r "cd('$xASLdir');ExploreASL_Master('$DataParPath', true, true, $i, $nParallel);system('screen -SX SABRE$i kill');" &
done