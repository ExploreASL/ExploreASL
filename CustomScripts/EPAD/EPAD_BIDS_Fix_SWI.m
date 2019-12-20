function EPAD_Fix_SWI_BIDS(AnalysisDir)
%EPAD_Fix_SWI_BIDS Fix dcm2nii SWI conversion per BIDS
%
% FORMAT: EPAD_Fix_SWI_BIDS(AnalysisDir)
% 
% INPUT:
%   AnalysisDir - path to folder containing the NIfTI/BIDS data (e.g. /data/RAD/share/EPAD/analysis)
%
% OUTPUT: n/a
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function manages the SWI BIDS structure, for all subjects.
%              It puts all SWI files in a single SWI folder & manages the
%              BIDS suffixes (i.e. run-1 run-2 etc, PartMag & PartPhase
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: EPAD_Fix_SWI_BIDS(AnalysisDir);
% __________________________________
% Copyright 2015-2019 ExploreASL

SubjectList = xASL_adm_GetFsList(AnalysisDir,'^\d{3}EPAD\d*$', true, [], [], [0 Inf]);

fprintf('%s','Adjusting SWI BIDS:  0%');

ScanType= 'swi';
DicomStr = 'DummyDicom_';

for iSubject=1:length(SubjectList)
    xASL_TrackProgress(iSubject,length(SubjectList));
    SWIdir = fullfile(AnalysisDir, SubjectList{iSubject}, 'swi_part_mag_1');
    
    Flist{1} = xASL_adm_GetFileList(SWIdir, '^swi_part_mag(|_run(-|_)\d*).*\.(nii|nii\.gz)$', 'FPList', [0 Inf]);
    Flist{2} = xASL_adm_GetFileList(SWIdir, '^swi_part_mag(|_run(-|_)\d*).*\.json$', 'FPList', [0 Inf]);
    Flist{3} = xASL_adm_GetFileList(SWIdir, '^swi_part_mag(|_run(-|_)\d*).*_parms\.mat$', 'FPList', [0 Inf]);
    Flist{4} = xASL_adm_GetFileList(SWIdir, ['^' DicomStr 'swi_part_mag(|_run(-|_)\d*).*\.dcm$'], 'FPList', [0 Inf]);
    
    % NB: _parms.mat only exists for the first nii
    
    for iFile=1:length(Flist)
        for iList=1:length(Flist{iFile})
            [DirOri, FileOri, ExtOri] = xASL_fileparts(Flist{iFile}{iList});

            % Get BIDS run index
            [Ind1, Ind2] = regexp(FileOri,'run(-|_)\d');
            RunString = FileOri(Ind1:Ind2);
            if isempty(RunString)
                RunString = 'run-1';
            end
            RunString(4) = '-';

            if iList==1
                ContrastString = 'PartMag';
            elseif iList==2
                ContrastString = 'PartPhase';
            end

            DirNew = fullfile(AnalysisDir, SubjectList{iSubject}, ScanType);
            DestName = [ScanType '_' RunString '_' ContrastString];

            PathOri = fullfile(DirOri, [FileOri ExtOri]);
            if strcmp(ExtOri,'.mat')
                ExtOri = '_parms.mat';
            elseif strcmp(ExtOri,'.dcm')
                DestName = [DicomStr DestName];
            end
            PathDest = fullfile(DirNew, [DestName ExtOri]);

            if exist(PathOri,'file')
                xASL_adm_CreateDir(DirNew);
                xASL_Move(PathOri, PathDest, true);
            end
        end
    end
    
    if exist(SWIdir,'dir') && isempty(xASL_adm_GetFileList(SWIdir, '.*', 'FPListRec', [0 Inf]))
       rmdir(SWIdir);
    end
end

fprintf('\n');

end