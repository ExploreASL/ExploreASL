#!/bin/bash
# Script to run all weekly tests in order
# Turn on/off the booleans to activate certain tests

# Set variables for Workerscript
MyScript=/scratch/mshammer/ExploreASL/Testing/xASL_test_RunInstance.sh
# set command to current matlab
Matlab=matlab-R2019b 

# Define folders.
scratchDir=/scratch/radv/mshammer
XASLDIR=${scratchDir}/ExploreASL
ReferenceTSV=${XASLDIR}/Testing/Reference/ReferenceValues.tsv
TestDataSetSourceDir=${scratchDir}/TestDataSets
TestDataSetWorkspaceDir=${TMP}/TestDataSetWorkspace
UnitTestingDir=${scratchDir}/Testing
FlavorTestConfig=${scratchDir}/WeeklyTest/FlavorConfig.json
ResultMasterDir=${scratchDir}/WeeklyTest/Results

# Email adress to send results to
EmailAdress=m.hammer@amsterdamumc.nl

# Booleans
bPull=false
bSPMTest=false
bUnitTest=true
bFlavorTest=true
bTestDataSet=true
bCompile=true
bSummary=true
bEmail=false
bParallel=false



# Make the results directory timed conform ISO 8601
today=$(date +"%FT%H:%M%:z") 
ResultDirToday=${ResultMasterDir}/${today}
mkdir ${ResultDirToday}

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
echo "We're testing version ${RepositoryVersion}"

# Run SPM test (no output?, fix that)
if ${bSPMTest}; then
	echo "bSPMTest was ${bSPMTest}"
	cd ${XASLDIR}
    nice -n 10 ${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_test_SPM('${TestDataSetWorkspaceDir}', false);exit;"
fi

# Run UnitTest
if ${bUnitTest}; then
	echo "bUnitTest was ${bUnitTest}"
	cd ${XASLDIR}
    nice -n 10 ${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_test_UnitTesting(false);exit;"
	mv ${UnitTestingDir}/*results.mat ${ResultDirToday}
	mv ${UnitTestingDir}/*comparison.tsv ${ResultDirToday}
fi

# Run Flavor Test Parallelize?
if ${bFlavorTest}; then
	echo "bFlavorTest was ${bFlavorTest} "
	cd ${XASLDIR}
    nice -n 10 ${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();config=spm_jsonread('${FlavorTestConfig}');xASL_test_Flavors(config, false, false);exit;"
	mv ${XASLDIR}/Testing/*results.mat ${ResultDirToday}
	mv ${XASLDIR}/Testing/*comparison.tsv ${ResultDirToday}
fi

if ${bTestDataSet}; then
	echo "bTestDataSet was ${bTestDataSet} "
	# Run Explore ASL on the TestDataSet Directory
	# Copy to a reference location and go there
	rm -r ${TestDataSetWorkspaceDir}
	cp -R ${TestDataSetSourceDir} ${TestDataSetWorkspaceDir} 
	cd ${TestDataSetWorkspaceDir}

	# create an array of all folders in the reference directory
	FolderArray=(*/)
	lengthDir="${#FolderArray[@]}"
	cd ${XASLDIR}

	# automatic status and lock file removal to ensure all steps are run
	echo "removing lock and status folders from ${TestDataSetWorkspaceDir}/*/derivatives/ExploreASL/lock/*/*/*/locked"
	for ((i=0; i<${lengthDir}; i++));
	do
		rm -fd ${TestDataSetWorkspaceDir}/${FolderArray[i]}derivatives/ExploreASL/lock/*/*/*/locked;
		rm -f ${TestDataSetWorkspaceDir}/${FolderArray[i]}derivatives/ExploreASL/lock/*/*/*/*.status 
	done

	# DEBUGGING
	# Set to true if you want an interactive slurm shell to test if the environmental variables are passed on correctly. 
	if false; then
		srun --export=ALL,WORKER=1,NWORKERS=${lengthDir},{XASLDIR}=${XASLDIR},DATAFOLDER=${TestDataSetWorkspaceDir}${FolderArray[0]} --job-name 'interactive-shell' --cpus-per-task 2 --mem-per-cpu 8 --time 0:00:30 --pty bash
	fi

	if bParallel; then
		# Make a seperate slurm job for every folder in the test directory to speed up process.
		for (( i=0; i<${lengthDir}; i++ ));
		do
			sbatch -W --export=ALL,WORKER=${i},NWORKERS=${lengthDir},XASLDIR=${XASLDIR},DATAFOLDER=${TestDataSetWorkspaceDir}${FolderArray[i]} --job-name 'xASL_'.${i}.${lengthDir}.run $MyScript
		done
		wait 
	# alternatively 
	# while `squeue -u mshammer | wc -l` >= 2; do 
	# sleep 5m
	# done
	else
		for (( i=0; i<${lengthDir}; i++ ));
		do
			nice -n 10 ${Matlab} -nodesktop -nosplash -r "cd('$XASLDIR');ExploreASL('${TestDataSetWorkspaceDir}/${FolderArray[i]}', 0, 1);exit;"
		done
	fi

	# Debugging, dont run this yet
	if true; then
		# Compare results to Reference Values
		nice -n 10 `${Matlab} -nodesktop -nosplash -r "cd('${XASLDIR}');ExploreASL();xASL_test_CompareReference('${ReferenceTSV}','${TestDataSetWorkspaceDir}');exit;"`
		mv ${TestDataSetWorkspaceDir}/*.tsv ${ResultDirToday}
		mv ${TestDataSetWorkspaceDir}/*ResultsTable.mat ${ResultDirToday}
	fi

	# Clean up temporary files
	rm -r ${TestDataSetWorkspaceDir}
fi

if ${bCompile}; then
	mkdir ${ResultDirToday}/compilation
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
		mail -s 'xASL git commit detected' -a ${RepositoryVersion}_Results.txt m.hammer@amsterdamumc.nl <<< 'Git commit ${RepositoryVersion} resulted in changes in the test results.\n Changes are attached in the text file.' 
		exit 0
	else 
		echo "No Changes detected, aborting"
		rm ${RepositoryVersion}_Results.txt 
		exit 0
	fi 
fi