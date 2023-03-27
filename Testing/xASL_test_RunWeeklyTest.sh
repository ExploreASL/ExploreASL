#!/bin/bash
# Script to run all weekly tests.
# Turn on/off the booleans to activate certain tests.

# set command to current matlab
Matlab=matlab-R2019b 

# Define folders.
scratchDir=/scratch/radv/mshammer
XASLDIR=${scratchDir}/ExploreASL
ReferenceTSV=${XASLDIR}/Testing/Reference/ReferenceValues.tsv
FlavorTestConfig=${scratchDir}/WeeklyTest/FlavorConfig.json
FlavorDir=${scratchDir}/FlavorDatabase
TestDataSetSourceDir=${scratchDir}/TestDataSets
UnitTestingDir=${scratchDir}/Testing
ResultMasterDir=${scratchDir}/WeeklyTest/Results

# Temporary Folders, ALL CONTENT WILL BE REMOVED FROM THIS FOLDER.
TestDataSetWorkspaceDir=${TMP}/TestDataSetWorkspace

# Email adress to send results to
EmailAdress=m.hammer@amsterdamumc.nl

# Booleans
bPull=false
bSPMTest=true
bUnitTest=true
bFlavorTest=true
bTestDataSet=true
bCompile=false
bSummary=true
bEmail=false
iNiceness=10 # 0 in testing, 10 at weekend.

# Make the results directory timed conform ISO 8601
today=$(date +"%FT%H:%M%:z") 
ResultDirToday=${ResultMasterDir}/${today}
mkdir ${ResultDirToday}
VersionFile=${ResultDirToday}/VersionsTested.txt
LogFile=${ResultDirToday}/LogFile.txt
touch ${VersionFile}

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
fi

# Run UnitTest
if ${bUnitTest}; then
	cd ${UnitTestingDir}
	git pull
	UnitVersion=`git rev-parse --short HEAD` 
	cd ${XASLDIR}
	echo "Unit test directory was tested on version ${UnitVersion}." >>  ${VersionFile}

    nice -n ${iNiceness} ${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_test_UnitTesting(false);exit;"
	mv ${UnitTestingDir}/*results.mat ${ResultDirToday}
	mv ${UnitTestingDir}/*comparison.tsv ${ResultDirToday}
fi

# Run Flavor Test Parallelize?
if ${bFlavorTest}; then
	cd ${FlavorDir}
	# TODO : remove rawdata and Derivatives folders for clean rerun.
	git pull
	FlavorVersion=`git rev-parse --short HEAD` 
	cd ${XASLDIR}
	echo "Flavor database test directory was tested  on version ${FlavorVersion}." >>  ${VersionFile}

    nice -n ${iNiceness} ${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_test_Flavors();exit;"
	mv ${XASLDIR}/Testing/*results.mat ${ResultDirToday}
	mv ${XASLDIR}/Testing/*comparison.tsv ${ResultDirToday}
	# TODO : remove rawdata and Derivatives folders for clean rerun.
fi

if ${bTestDataSet}; then
	# Run Explore ASL on the TestDataSet Directory
	# Copy to a reference location and go there
	cd ${TestDataSetSourceDir}
	git pull
	TestDataSet=`git rev-parse --short HEAD` 
	cd ${XASLDIR}
	rm -rf ${TestDataSetWorkspaceDir}
	cp -R ${TestDataSetSourceDir} ${TestDataSetWorkspaceDir} 
	cd ${TestDataSetWorkspaceDir}
	echo "TestDataSet directory was tested on version ${bTestDataSet}." >> ${VersionFile}

	# create an array of all folders in the reference directory
	FolderArray=(*/)
	lengthDir="${#FolderArray[@]}"
	cd ${XASLDIR}

	# Run all test
	for (( i=0; i<${lengthDir}; i++ ));
	do
		nice -n ${iNiceness} ${Matlab} -nodesktop -nosplash -r "cd('$XASLDIR');ExploreASL();ExploreASL('${TestDataSetWorkspaceDir}/${FolderArray[i]}', 0, 1, false);exit;"
	done

	# Compare results to Reference Values
	nice -n ${iNiceness} ${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_test_CompareReference('${ReferenceTSV}','${TestDataSetWorkspaceDir}');exit;"
	mv ${TestDataSetWorkspaceDir}/*.tsv ${ResultDirToday}
	mv ${TestDataSetWorkspaceDir}/*ResultsTable.mat ${ResultDirToday}

	# Clean up temporary files
	rm -rf ${TestDataSetWorkspaceDir}
fi

if ${bCompile}; then
	mkdir ${ResultDirToday}/compilation
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
		mail -s 'xASL git commit detected' -a ${RepositoryVersion}_Results.txt m.hammer@amsterdamumc.nl <<< 'Git commit ${RepositoryVersion} resulted in changes in the test results.\n Changes are attached in the text file.' 
		exit 0
	else 
		echo "No Changes detected, aborting"
		rm ${RepositoryVersion}_Results.txt 
		exit 0
	fi 
fi
