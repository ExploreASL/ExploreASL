function xASL_qc_ObtainQCCategoriesFromJPG(x)
%xASL_qc_ObtainQCCategoriesFromJPG Obtain QC categories as covariant/set
%
% FORMAT: xASL_qc_ObtainQCCategoriesFromJPG(x)
% 
% INPUT:
%   x                            - struct containing statistical pipeline environment parameters (REQUIRED)
%
% OUTPUT: n/a
% -------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function obtains QC categories as covariant/set,
%              based on the JPGs in //Population/ASLCheck. These are initially sorted by
%              spatial CoV, and should be visually checked & put in the correct folder.
% -------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_qc_ObtainQCCategoriesFromJPG(x);
% __________________________________
% Copyright 2015-2019 ExploreASL


%% Find QC folders
QCnames = {'1_CBFContrast' '2_VascularContrast' '3_ArtifactContrast' '4_Unknown_sCoV'};
iNext = 1;
for iQC=1:length(QCnames)
    DirName{iQC} = fullfile(x.D.ASLCheckDir, QCnames{iQC});
    if exist(DirName{iQC},'dir')
        FileList = xASL_adm_GetFileList(DirName{iQC}, 'Tra_qCBF.*jpg$', 'List',[0 Inf]);
        for iList=1:length(FileList)
            [StartInd, EndInd] = regexp(FileList{iList},'ASL_\d+');
            SessionID = FileList{iList}(StartInd:EndInd);
            SubjectID = FileList{iList}(length('Tra_qCBF_')+1:StartInd-2);
            QCcategory{iNext,1} = SubjectID;
            QCcategory{iNext,2} = SessionID;
            QCcategory{iNext,3} = QCnames{iQC};
            iNext = iNext+1;
        end
    elseif iQC<4
        warning(['Category folder missing: ' DirName{iQC}]);
    end
end
            
QC_SavePath = fullfile(x.dir.xASLDerivatives, 'QCcategory.mat');
save(QC_SavePath,'QCcategory');

fprintf('\n');

end
    