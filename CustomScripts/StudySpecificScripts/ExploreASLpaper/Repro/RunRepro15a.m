function RunRepro15a
%RunRepro Summary of this function goes here
%   Detailed explanation goes here

addpath('c:\ExploreASL\CustomScripts\ExploreASLpaper');
cd('c:\ExploreASL');

ExploreASL_Master('C:\BackupWork\ASL\ReproWork\EPAD_Win15a_New_1\DataPar.m',true,true);
ExploreASL_Master('C:\BackupWork\ASL\ReproWork\Sleep_Win15a_New_1\DataPar.m',true,true);
ExploreASL_Master('C:\BackupWork\ASL\ReproWork\Novice_Win15a_New_1\DataPar.m',true,true);

end

