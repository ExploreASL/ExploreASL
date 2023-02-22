#!/bin/bash
# Script to run all weekly tests.
# Turn on/off the booleans to activate certain tests.

<<<<<<< HEAD
<<<<<<< HEAD
=======
# Set variables for Workerscript
MyScript=/scratch/mshammer/ExploreASL/Testing/xASL_test_RunInstance.sh
>>>>>>> 4690c8f9 ( #1247 added Version summary and pulling/resetting flavor and dataset to clean)
=======
>>>>>>> f6c0e0ae ( #1247 Removed paralels for now)
# set command to current matlab
Matlab=matlab-R2019b 

# Define folders.
scratchDir=/scratch/radv/mshammer
<<<<<<< HEAD
XASLDIR=${scratchDir}/WeeklyTest/ExploreASL
ReferenceTSV=${XASLDIR}/Testing/Reference/ReferenceValues.tsv
FlavorTestConfig=${scratchDir}/WeeklyTest/FlavorConfig.json
FlavorDir=${scratchDir}/WeeklyTest/FlavorDatabase
TestDataSetSourceDir=${scratchDir}/WeeklyTest/TestDataSets
UnitTestingDir=${scratchDir}/WeeklyTest/Testing
ResultMasterDir=${scratchDir}/WeeklyTest/Results

# Temporary Folders, ALL CONTENT WILL BE REMOVED FROM THIS FOLDER.
=======
XASLDIR=${scratchDir}/ExploreASL
ReferenceTSV=${XASLDIR}/Testing/Reference/ReferenceValues.tsv
FlavorTestConfig=${scratchDir}/WeeklyTest/FlavorConfig.json
FlavorDir=${scratchDir}/FlavorDatabase
TestDataSetSourceDir=${scratchDir}/TestDataSets
UnitTestingDir=${scratchDir}/Testing
ResultMasterDir=${scratchDir}/WeeklyTest/Results

# Temporary Folders.
>>>>>>> 4690c8f9 ( #1247 added Version summary and pulling/resetting flavor and dataset to clean)
TestDataSetWorkspaceDir=${TMP}/TestDataSetWorkspace

# Email adress to send results to
EmailAdress=m.hammer@amsterdamumc.nl

# Booleans
bPull=false
bSPMTest=true
bUnitTest=true
bFlavorTest=true
bTestDataSet=true
<<<<<<< HEAD
bCompile=false
bSummary=true
bEmail=false
iNiceness=10 # 0 in testing, 10 at weekend.
=======
bCompile=true
bSummary=true
bEmail=false
<<<<<<< HEAD
bParallel=false
>>>>>>> 4690c8f9 ( #1247 added Version summary and pulling/resetting flavor and dataset to clean)
=======
iNiceness=10
>>>>>>> f6c0e0ae ( #1247 Removed paralels for now)

# Make the results directory timed conform ISO 8601
today=$(date +"%FT%H:%M%:z") 
ResultDirToday=${ResultMasterDir}/${today}
mkdir ${ResultDirToday}
VersionFile=${ResultDirToday}/VersionsTested.txt
<<<<<<< HEAD
LogFile=${ResultDirToday}/LogFile.txt
touch ${VersionFile}
=======
touch VersionFile
>>>>>>> 4690c8f9 ( #1247 added Version summary and pulling/resetting flavor and dataset to clean)

# Initialize some variables
cd ${XASLDIR}

# Fetch remote branch, see if up to date and pull changes.
if ${bPull}; then
	echo "bPull was ${bPull}"
	cd ${XASLDIR}
	git remote update
	if git status -uno | grep -q 'Your branch is up-to-date'; then
		echo "Branch up to date, quitting CDCI" 
		exit 0
	elif git status -uno | grep -q 'Your branch is behind'; then
		echo "Branch delayed, Proceeding to CDCI"
		git pull
	else 
		echo "Error, Branch status not found, quitting CDCI."
		exit 1
	fi
fi

# Get current Commit hash
RepositoryVersion=`git rev-parse --short HEAD` 
<<<<<<< HEAD
echo "We're testing version on ExploreASL version ${RepositoryVersion}." >>  ${VersionFile}

# Run SPM test (no output?, fix that)
if ${bSPMTest}; then
	# Run Explore ASL on the TestDataSet Directory
	# Copy to a reference location and go there
	cd ${TestDataSetSourceDir}
	git pull
	rm -rf ${TestDataSetWorkspaceDir}
	cp -R ${TestDataSetSourceDir} ${TestDataSetWorkspaceDir} 
	echo "SPMTest was conducted in ExploreASL version ${RepositoryVersion}." >>  ${VersionFile}

	cd ${XASLDIR}
    nice -n ${iNiceness} ${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_test_SPM('${TestDataSetWorkspaceDir}', false);exit;"

	# Clean up temporary files
	rm -rf ${TestDataSetWorkspaceDir}
=======
echo "We're testing version on ExploreASL version ${RepositoryVersion}." >> VersionFile

# Run SPM test (no output?, fix that)
if ${bSPMTest}; then
	echo "bSPMTest was conducted in ExploreASL version ${RepositoryVersion}." >> VersionFile
	cd ${XASLDIR}
<<<<<<< HEAD
    nice -n 10 ${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_test_SPM('${TestDataSetWorkspaceDir}', false);exit;"
>>>>>>> 4690c8f9 ( #1247 added Version summary and pulling/resetting flavor and dataset to clean)
=======
    nice -n ${iNiceness} ${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_test_SPM('${TestDataSetWorkspaceDir}', false);exit;"
>>>>>>> f6c0e0ae ( #1247 Removed paralels for now)
fi

# Run UnitTest
if ${bUnitTest}; then
	cd ${UnitTestingDir}
	git pull
	UnitVersion=`git rev-parse --short HEAD` 
	cd ${XASLDIR}
<<<<<<< HEAD
	echo "Unit test directory was tested on version ${UnitVersion}." >>  ${VersionFile}

    nice -n ${iNiceness} ${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_test_UnitTesting(false);exit;"
=======
	echo "Unit test directory was tested on version ${UnitVersion}." >> VersionFile
<<<<<<< HEAD
    nice -n 10 ${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_test_UnitTesting(false);exit;"
>>>>>>> 4690c8f9 ( #1247 added Version summary and pulling/resetting flavor and dataset to clean)
=======
    nice -n ${iNiceness} ${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_test_UnitTesting(false);exit;"
>>>>>>> f6c0e0ae ( #1247 Removed paralels for now)
	mv ${UnitTestingDir}/*results.mat ${ResultDirToday}
	mv ${UnitTestingDir}/*comparison.tsv ${ResultDirToday}
fi

# Run Flavor Test Parallelize?
if ${bFlavorTest}; then
	cd ${FlavorDir}
	git pull
	FlavorVersion=`git rev-parse --short HEAD` 
	cd ${XASLDIR}
<<<<<<< HEAD
	echo "Flavor database test directory was tested  on version ${FlavorVersion}." >>  ${VersionFile}

    nice -n ${iNiceness} ${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_test_Flavors();exit;"
	mv ${XASLDIR}/Testing/*results.mat ${ResultDirToday}
	mv ${XASLDIR}/Testing/*comparison.tsv ${ResultDirToday}
	mv ${XASLDIR}/Testing/*flavor_loggingtable.tsv ${ResultDirToday}
=======
	echo "Flavor database test directory was tested  on version ${FlavorVersion}." >> VersionFile
    nice -n ${iNiceness} ${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();config=spm_jsonread('${FlavorTestConfig}');xASL_test_Flavors(config, false, false);exit;"
	mv ${XASLDIR}/Testing/*results.mat ${ResultDirToday}
	mv ${XASLDIR}/Testing/*comparison.tsv ${ResultDirToday}
>>>>>>> 4690c8f9 ( #1247 added Version summary and pulling/resetting flavor and dataset to clean)
fi

if ${bTestDataSet}; then
	# Run Explore ASL on the TestDataSet Directory
	# Copy to a reference location and go there
	cd ${TestDataSetSourceDir}
	git pull
	TestDataSet=`git rev-parse --short HEAD` 
	cd ${XASLDIR}
<<<<<<< HEAD
	rm -rf ${TestDataSetWorkspaceDir}
	cp -R ${TestDataSetSourceDir} ${TestDataSetWorkspaceDir} 
	cd ${TestDataSetWorkspaceDir}
	echo "TestDataSet directory was tested on version ${bTestDataSet}." >> ${VersionFile}
=======
	echo "TestDataSet directory was tested on version ${bTestDataSet}." >> VersionFile
	rm -rf ${TestDataSetWorkspaceDir}
	cp -R ${TestDataSetSourceDir} ${TestDataSetWorkspaceDir} 
	cd ${TestDataSetWorkspaceDir}
>>>>>>> 4690c8f9 ( #1247 added Version summary and pulling/resetting flavor and dataset to clean)

	# create an array of all folders in the reference directory
	FolderArray=(*/)
	lengthDir="${#FolderArray[@]}"
	cd ${XASLDIR}

<<<<<<< HEAD
	# Run all test
	for (( i=0; i<${lengthDir}; i++ ));
	do
		nice -n ${iNiceness} ${Matlab} -nodesktop -nosplash -r "cd('$XASLDIR');ExploreASL();ExploreASL('${TestDataSetWorkspaceDir}/${FolderArray[i]}', 0, 1, false);exit;"
	done

	# Compare results to Reference Values
	nice -n ${iNiceness} ${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_test_CompareReference('${ReferenceTSV}','${TestDataSetWorkspaceDir}');exit;"
	mv ${TestDataSetWorkspaceDir}/*.tsv ${ResultDirToday}
	mv ${TestDataSetWorkspaceDir}/*ResultsTable.mat ${ResultDirToday}
=======
	# automatic status and lock file removal to ensure all steps are run
	echo "removing lock and status folders from ${TestDataSetWorkspaceDir}/*/derivatives/ExploreASL/lock/*/*/*/locked"
	for ((i=0; i<${lengthDir}; i++));
	do
		rm -fd ${TestDataSetWorkspaceDir}/${FolderArray[i]}derivatives/ExploreASL/lock/*/*/*/locked;
		rm -f ${TestDataSetWorkspaceDir}/${FolderArray[i]}derivatives/ExploreASL/lock/*/*/*/*.status 
	done

	# Run all test
	for (( i=0; i<${lengthDir}; i++ ));
	do
		nice -n ${iNiceness} ${Matlab} -nodesktop -nosplash -r "cd('$XASLDIR');ExploreASL('${TestDataSetWorkspaceDir}/${FolderArray[i]}', 0, 1);exit;"
	done


	# Debugging, dont run this yet
	if true; then
		# Compare results to Reference Values
		nice -n 10 `${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_test_CompareReference('${ReferenceTSV}','${TestDataSetWorkspaceDir}');exit;"`
		mv ${TestDataSetWorkspaceDir}/*.tsv ${ResultDirToday}
		mv ${TestDataSetWorkspaceDir}/*ResultsTable.mat ${ResultDirToday}
	fi
>>>>>>> 4690c8f9 ( #1247 added Version summary and pulling/resetting flavor and dataset to clean)

	# Clean up temporary files
	rm -rf ${TestDataSetWorkspaceDir}
fi

if ${bCompile}; then
	mkdir ${ResultDirToday}/compilation
<<<<<<< HEAD
	nice -n ${iNiceness} ${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_dev_MakeStandalone('${ResultDirToday}/compilation','${TestDataSetWorkspaceDir}');exit;"
fi

if ${bSummary}; then
	nice -n ${iNiceness} ${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_test_Summarize('${ResultDirToday}');exit;"
fi

if ${bEmail}; then
	echo "bEmail was ${bEmail}" >> ${VersionFile}
	cd ${ResultDirToday}

	# Email when outputfile exists and is not empty, delete empty output if exists
	if [ -s ${RepositoryVersion}_Results.txt ]; then
		echo "Sending email to ${EmailAdress}" >> ${VersionFile}
=======
	nice -n 10 `${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_dev_MakeStandalone('${ResultDirToday}/compilation','${TestDataSetWorkspaceDir}');exit;"`
fi

if ${bSummary}; then
	nice -n 10 `${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_test_Summarize('${ResultDirToday}');exit;"`
fi

if ${bEmail}; then
	echo "bEmail was ${bEmail} "
	cd ${ResultDirToday}
	# Email when outputfile exists and is not empty, delete empty output if exists
	if [ -s ${RepositoryVersion}_Results.txt ]; then
		echo "Sending email to ${EmailAdress}"
>>>>>>> 4690c8f9 ( #1247 added Version summary and pulling/resetting flavor and dataset to clean)
		mail -s 'xASL git commit detected' -a ${RepositoryVersion}_Results.txt m.hammer@amsterdamumc.nl <<< 'Git commit ${RepositoryVersion} resulted in changes in the test results.\n Changes are attached in the text file.' 
		exit 0
	else 
		echo "No Changes detected, aborting"
		rm ${RepositoryVersion}_Results.txt 
		exit 0
	fi 
<<<<<<< HEAD
fi
=======
fi
>>>>>>> 4690c8f9 ( #1247 added Version summary and pulling/resetting flavor and dataset to clean)
