# note to first remove any existing locked folders
let nParallel=10;
MatlabPath=/usr/local/apps/matlab/R2018b/bin/matlab;
DataFolder=/radshare/SABRE/Analysis2; 
DataParPath=$DataFolder/DataParameters_HiQ.json;
rm $DataFolder/lock/*/*/*/locked -d; # remove locked folders
# srrm -rf $DataFolder/lock/xASL_module_Population; # Population lock folder
xASLdir=/radshare/ExploreASL_Test/ExploreASL;
cd $xASLdir;
for (( i=1; i<=$nParallel; i++ ));
do screen -dmSL SABRE$i nice -n 10 $MatlabPath -nodesktop -nosplash -r "cd('$xASLdir');ExploreASL_Master('$DataParPath', true, true, $i, $nParallel);system('screen -SX SABRE$i kill');" &
done


# note to remove locks