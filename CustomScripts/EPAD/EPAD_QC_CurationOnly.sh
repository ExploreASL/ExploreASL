MatlabPath=$(which matlab); # find Matlab path
ROOT=/mnt/s4e_data/RAD/share/EPAD/EPAD500.1; # dir of data (should have /raw inside, or will be created)
ASLdir=/mnt/s4e_data/RAD/share/EPAD/scripts/ExploreASL; # dir of ExploreASL
cd $xASLdir
cd CustomScripts/EPAD

screen -dmS CurationOnly nice -n 10 $MatlabPath -nodesktop -nosplash -r "addpath(pwd);cd ../..;EPAD_QC_Master('$ROOT',0,0,0)"
# this code is explained in parts below:
# screen -dmS CurationOnly -> creates a screen with this name, to run this in the background
# nice -n 10 -> avoid too high priority
# $MatlabPath -nodeskopt -nosplash -r -> start up Matlab in CLI mode without start up screen, running matlab code after -r 
# "addpath(pwd);cd ../..;" -> this manages paths
# EPAD_QC_Master -> this contains the full QC pipeline, but will skip parts that have been run before
# 0,0 -> don't use DCMTK, don't check permissions
# last 0 means -> don't run ExploreASL (last step) -> so this will only run the curation & dcm2niiX, not ExploreASL image processing & QC