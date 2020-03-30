let nParallel=6;

MatlabPath=/Applications/MATLAB_R2016b.app/bin/matlab; # find Matlab path
DataFolder=/Users/henk/ExploreASL/ASL/Test_Twins2Run/LST_CleanUp_Yes_HJ_CleanUp_No;
DataParPath=$DataFolder/DATA_PAR.json;
xASLdir=/Users/henk/ExploreASL/ExploreASL;
cd $xASLdir

rm $DataFolder/lock/*/*/*/locked -d; # remove locked folders
rm -rf $DataFolder/lock/xASL_module_Population; # Population lock folder

for (( i=1; i<=$nParallel; i++ ));
do screen -dmSL TwinsTrial$i nice -n 10 $MatlabPath -nodesktop -nosplash -r "cd('$xASLdir');ExploreASL_Master('$DataParPath', true, true, $i, $nParallel);exit;system('screen -SX TwinsTrial$i kill');" &
done