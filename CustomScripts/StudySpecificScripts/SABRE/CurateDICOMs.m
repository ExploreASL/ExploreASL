%% Curate SABRE data
% Obtain regular expression

ExploreASL_Master('',0);

AnalysisDir = '/Users/henk/surfdrive/SABRE/analysis';
BaMoSDir = '/Users/henk/surfdrive/SABRE/Carole/SABRE';

RegExp = '\d{5}\d*';

xASL_adm_ImportBaMoS(AnalysisDir, BaMoSDir, 0, RegExp);


    % for EPAD

    % here we check if the scanner prefix is a site prefix (e.g. 010 instead of 110)
    % & also if there is a '-' instead of 'EPAD'
    Name1 = ['0' NameNew(2:3) '-' NameNew(8:end)];
    Name2 = NameNew(8:end);
    
    
    
    % 501609 -> two different persons? -> check with who this one belongs
    % 191658 -> second was incorrectly rotated, removed
    % 133575 -> 2nd smaller FoV, removed

    
    
% 1) Wait for DICOM curation to be ready
% 2) Run dcm2niiX
% 3) Check 501609 -> two different persons
% 4) Rerun WMH cloning
% 5) Check which ones have ASL, exclude those without
% 6) Start processing structural data
% 7) Fix DataPar ASL parameters
% 8) Jantine -> Visual QC
% 9) Jantine -> correct for Hct (but don't use immediately, check individual vs group)