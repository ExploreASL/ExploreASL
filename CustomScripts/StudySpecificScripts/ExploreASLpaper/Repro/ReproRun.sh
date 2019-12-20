# DirExploreASL=/mnt/c/Backup/OtherToolboxesEtc/ExploreASL_Repro/ExploreASL
DirExploreASL=/mnt/c/ExploreASL
# DirMatlab18b=/mnt/c/Program\ Files/MATLAB/R2018a/bin
# DirMatlab14a=/mnt/c/Program\ Files/MATLAB/R2014a/bin
# DirMatlab11a=/mnt/c/Program\ Files/MATLAB/R2011a/bin

if [ ! -d $DirExploreASL ]
then
	echo "warning: ExploreASL folder did not exist:" $DirExploreASL
	exit
fi
	

cd $DirExploreASL

# Step 1: Create TestDataSet_Processed if it doesn't exist
if [ -d "./TestDataSet_Processed" ]
then
	echo "TestDataSet_Processed already existed"
	
	# Step 2: Test other versions
	echo "Make sure to set DATA_PAR.m in TestDataSet to x.Quality==1 & bRepro==1"

	echo "Copying & running other datasets"
	cp -rf ./TestDataSet ./TestDataSet_Win2014a
	cp -rf ./TestDataSet ./TestDataSet_Win2011a

	/mnt/c/Program\ Files/MATLAB/R2014a/bin/matlab.exe -nodesktop -nosplash -r "ExploreASL_Master('c:\ExploreASL\TestDataSet_Win2014a\DATA_PAR.m', true, true);exit"
	/mnt/c/Program\ Files/MATLAB/R2011a/bin/matlab.exe -nodesktop -nosplash -r "ExploreASL_Master('c:\ExploreASL\TestDataSet_Win2011a\DATA_PAR.m', true, true);exit"	
	
else
	echo "Creating & running TestDataSet_Processed"
	cp -rf ./TestDataSet ./TestDataSet_Processed
	/mnt/c/Program\ Files/MATLAB/R2018b/bin/matlab.exe -nodesktop -nosplash -r "ExploreASL_Master('c:\ExploreASL\TestDataSet_Processed\DATA_PAR.m', true, true);exit"
fi


