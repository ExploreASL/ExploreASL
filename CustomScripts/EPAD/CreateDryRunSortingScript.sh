# source=c:\Backup\ASL\EPAD\raw2
# dest=c:\Backup\ASL\EPAD\test
# ExploreASLdir=c:\ExploreASL

# [ -d $dest ] || mkdir $dest # create destination root folder
# rsync -av -f"+ */" -f"- *" "$source" "$dest" # copy folder structure

# Run matlab to copy 1 file per folder
/mnt/c/'Program Files'/MATLAB/R2018b/bin/matlab.exe -nojvm -nodesktop -nosplash -r "CreateDryRunSortingScript;"
