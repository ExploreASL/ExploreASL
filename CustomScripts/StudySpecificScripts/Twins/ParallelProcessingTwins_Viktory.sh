let nParallel=10;
MatlabPath=/usr/local/apps/matlab/R2018b/bin/matlab;
DataFolder=/radshare/Twins/Twins_ExploreASL; 
DataParPath=$DataFolder/DATA_PAR_TwinStudy.json;
rm $DataFolder/lock/*/*/*/locked -d; # remove locked folders
rm -rf $DataFolder/lock/xASL_module_Population; # Population lock folder
xASLdir=/radshare/EPAD/scripts/ExploreASL;
cd $xASLdir
for (( i=1; i<=$nParallel; i++ ));
do screen -dmSL TWINS$i nice -n 10 $MatlabPath -nodesktop -nosplash -r "cd('$xASLdir');ExploreASL_Master('$DataParPath', true, true, $i, $nParallel);system('screen -SX TWINS$i kill');" &
done