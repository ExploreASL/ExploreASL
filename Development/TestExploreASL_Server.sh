# set folders
xASLdir=/radshare/ExploreASL_Test/ExploreASL;
MatlabPath=/usr/local/apps/matlab/R2018b/bin/matlab;
TestDataPath=/radshare/ExploreASL_Test/ExploreASL_TestCases;
ProcessedDataPath=/radshare/ExploreASL_Test/ExploreASL_TestCasesProcessed;
cd $xASLdir;
# get most recent ExploreASL edits from git
git fetch --all;
git reset --hard origin/master;
# run matlab
$MatlabPath -nodesktop -nosplash -r "cd('$xASLdir/Development');xASL_qc_TestExploreASL('$TestDataPath','$ProcessedDataPath',2,0,'$MatlabPath');"