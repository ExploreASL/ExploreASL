function xASL_adm_ImportBaMoS_NewStructure(AnalysisDir, BaMoSDir, bPullPush, RegExp)
%xASL_adm_Import_BaMoS_NewStructure Prepare ASL parameter files for different EPAD vendors/sites
%
% FORMAT: xASL_adm_ImportBaMoS(AnalysisDir, BaMoSDir, 0, RegExp)
% 
% INPUT:
%   ROOT        - path to root folder containing the /analysis folder with NIfTI/BIDS data (REQUIRED)
%   CaroleDir   - path to root folder of Carole' segmentations & FLAIRs
%                 (OPTIONAL, trying a previous location if not provided)
%   bPullPush   - whether to create the subject list from the AnalysisDir
%                 (false, pull) or from the BaMoSdir (true, push) (OPTIONAL,
%                 DEFAULT=false)
%   RegExp      - regular expression for the SubjectID (REQUIRED)
%
% OUTPUT: n/a
% OUTPUT FILE:
%    /ROOT/analysis/NotCopiedList.json - contains list of missing WMH
%                                        segmentations (where we couldnt find both the FLAIR/WMH_SEGM from
%                                        Carole)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function searches for Carole' FLAIRs & WMH_SEGM & copies them into the ExploreASL/BIDS compliant
%              directory structure, for further use in EPAD analyses & QC.
%              We copy the FLAIR in addition to the WMH, as the WMH
%              segmentation will be aligned with Carole' FLAIR but not
%              necessarily with ours. For ExploreASL, these files are
%              copied to:
%              /ROOT/analysis/SubjectName/FLAIR.nii[.gz]
%              /ROOT/analysis/SubjectName/WMH_SEGM.nii[.gz]
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE for MacOS: xASL_adm_ImportBaMoS('/Users/henk/ExploreASL/ASL/SABRE/analysis', '/Users/henk/ExploreASL/ASL/SABRE/Carole', false, '\d{5}\d*');
% EXAMPLE for server: xASL_adm_ImportBaMoS('/radshare/SABRE/analysis', '/radshare/SABRE/Carole', false, '\d{5}\d*');
% __________________________________
% Copyright 2015-2019 ExploreASL


if nargin<3 || isempty(bPullPush)
    bPullPush = 0;
end

Dir2 = fullfile(BaMoSDir, 'ResultsProcessing');
if exist(Dir2, 'dir')
     BaMoSDir = Dir2;
end

% Adjust the regular expression
if strcmp(RegExp(1),'^')
    RegExp = RegExp(2:end);
end
if strcmp(RegExp(end),'$')
    RegExp = RegExp(1:end-1);
end

if bPullPush==0 % pull
    % use list from AnalysisDir, assuming that this contains the correct
    % subject names described by RegExp
    SubjList = xASL_adm_GetFileList(AnalysisDir, RegExp, 'List', [0 Inf], true);
else % push
    % use list from BaMosDir, doesn't have to have correct names
    SubjList = xASL_adm_GetFileList(BaMoSDir, ['.*' RegExp '.*\.nii'], 'FPListRec');
    for iSubj=1:length(SubjList)
        [StartInd, EndInd] = regexp(SubjList{iSubj}, RegExp);
        SubjList{iSubj} = SubjList{iSubj}(StartInd:EndInd);
    end
end
    
if isempty(SubjList)
    error('No subjects found');
end

NoFLAIRList = struct;
TooManyFLAIRList = struct;
NoWMHList = struct;
TooManyWMHList = struct;

fprintf('Importing FLAIR & WMH segmentations BaMoS:   ');
for iSubj=1:length(SubjList)
    xASL_TrackProgress(iSubj, length(SubjList));
    
    FLAIRnew = fullfile(AnalysisDir, SubjList{iSubj}, 'FLAIR.nii');
    WMHnew = fullfile(AnalysisDir, SubjList{iSubj}, 'WMH_SEGM.nii');
    
    FLAIRold = xASL_adm_GetFileList(BaMoSDir, ['FLAIR_' SubjList{iSubj} '\.nii'],'FPListRec', [0 Inf]);
    WMHold = xASL_adm_GetFileList(BaMoSDir, ['Lesion_' SubjList{iSubj} '\.nii'],'FPListRec', [0 Inf]);
    
    ExistFLAIR = false;
    ExistWMH = false;
    if length(FLAIRold)==0
        NoFLAIRList.(['n' num2str(length(fields(NoFLAIRList)))]) = SubjList{iSubj};
    elseif length(FLAIRold)>1
        TooManyFLAIRList.(['n' num2str(length(fields(TooManyFLAIRList)))]) = SubjList{iSubj};
        ExistFLAIR = true;
    else
        ExistFLAIR = true;
    end
    
    ExistWMH = false;
    if length(WMHold)==0
        NoWMHList.(['n' num2str(length(fields(NoWMHList)))]) = SubjList{iSubj};
    elseif length(FLAIRold)>1
        TooManyWMHList.(['n' num2str(length(fields(TooManyWMHList)))]) = SubjList{iSubj};
        ExistWMH = true;
    else
        ExistWMH = true;
    end    
    
    % The latest segmentation/FLAIR is usually the best
    if ExistFLAIR && ExistWMH
        xASL_Copy(FLAIRold{end}, FLAIRnew, true);
        xASL_Copy(WMHold{end}, WMHnew, true);
    end
end

fprintf('\n');

%% SAVE THE MISSING LIST HERE IN JSON FILES
ListsAre = {'NoFLAIRList' 'TooManyFLAIRList' 'NoWMHList' 'TooManyWMHList'};
for iList=1:length(ListsAre)
    if isstruct(eval(ListsAre{iList})) && ~isempty(fields(eval(ListsAre{iList})))
        SavePath = fullfile(AnalysisDir, [ListsAre{iList} '.json']);
        xASL_delete(SavePath);
        xASL_adm_SaveJSON(eval(ListsAre{iList}), SavePath);
    end

end


end