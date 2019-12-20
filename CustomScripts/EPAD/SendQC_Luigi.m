DstDir = '/scratch/hjmutsaerts/ForLuigi';
OriDir = '/scratch/hjmutsaerts/Data';
xASL_adm_CreateDir(DstDir);

xASL_adm_CopyMoveFileList(OriDir, DstDir, '^QC_collection.*\.json$', false, false, true, false, true);
cd(DstDir);
cd ..;
system('rsync -avz ForLuigi h.mutsaerts@145.121.126.53:/mnt/s4e_data/RAD/share/EPAD500_new/ForLuigi')