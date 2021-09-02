function xASL_imp_NII2BIDS_RunAnat(imPar, bidsPar, studyPar, subjectSessionLabel, outSessionPath, listRuns, iRun, nameSubjectSession)
%xASL_imp_NII2BIDS_RunAnat NII2BIDS conversion for a single sessions, single run.
%
% FORMAT: xASL_imp_NII2BIDS_RunAnat(bidsPar, studyPar, subjectSessionLabel, outSessionPath, listRuns, iRun)
% 
% INPUT:
% imPar               - JSON file with structure with import parameter (STRUCT, REQUIRED)
% bidsPar             - Output of xASL_imp_Config (STRUCT, REQUIRED)
% studyPar            - JSON file with the BIDS parameters relevant for the whole study (STRUCT, REQUIRED)
% subjectSessionLabel - subject-session label (CHAR ARRAY, REQUIRED)
% outSessionPath      - output session path (CHAR ARRAY, PATH, REQUIRED)
% listRuns            - list of runs (INTEGER, REQUIRED)
% iRun                - run number (INTEGER, REQUIRED)
% nameSubjectSession  - Name of subject & session (REQUIRED)
%
% OUTPUT:
% n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: NII2BIDS conversion for a single sessions, single run.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     n/a
% __________________________________
% Copyright 2015-2021 ExploreASL
    

    %% 1. Define the pathnames
    if ~isempty(listRuns)
        anatLabel = num2str(iRun);
        anatOutDirLabel = fullfile(outSessionPath,'anat');
        anatOutLabel = fullfile(outSessionPath,'anat',['sub-' subjectSessionLabel '_run-' num2str(iRun)]);
        anatOutLabelRelative = fullfile('anat',['sub-' subjectSessionLabel '_run-' num2str(iRun)]);
    else
        anatLabel = '';
        anatOutDirLabel = fullfile(outSessionPath,'anat');
        anatOutLabel = fullfile(outSessionPath,'anat',['sub-' subjectSessionLabel]);
        anatOutLabelRelative = fullfile('anat',['sub-' subjectSessionLabel]);
    end
    

    %% Iterate over files
    for iAnatType = bidsPar.listAnatTypes

        % Check if it exists
        anatPath = '';
		if xASL_exist(fullfile(imPar.TempRoot,nameSubjectSession,[iAnatType{1} '_1'],[iAnatType{1},'.nii']),'file')
			anatPath = fullfile(imPar.TempRoot,nameSubjectSession,[iAnatType{1} '_1'],iAnatType{1});
		elseif xASL_exist(fullfile(imPar.TempRoot,nameSubjectSession,[iAnatType{1},'.nii']),'file')
			anatPath = fullfile(imPar.TempRoot,nameSubjectSession,iAnatType{1});
		end

        % If anatomical file of this type exist, then BIDSify its structure
        if ~isempty(anatPath)

            % Create the anatomical directory
            xASL_adm_CreateDir(anatOutDirLabel);

            % Current scan name
            [~, scanName] = xASL_fileparts([anatOutLabel,'_' iAnatType{1}]);
            
            % Print anat name
            fprintf('Scan %s ...\n', scanName);
            
            % Move the NiFTI file
            anatNiiPath = [anatOutLabel,'_' iAnatType{1} '.nii.gz'];
            xASL_Move([anatPath '.nii'],anatNiiPath,1);

            % Load the JSON
            jsonAnat = spm_jsonread([anatPath,'.json']);

            % Save the JSON
            jsonAnat = xASL_bids_BIDSifyAnatJSON(jsonAnat,studyPar);
            jsonAnat = xASL_bids_VendorFieldCheck(jsonAnat);
            jsonAnat = xASL_bids_JsonCheck(jsonAnat,'');
            
            jsonWritePath = [anatOutLabel,'_' iAnatType{1} '.json'];
            spm_jsonwrite(jsonWritePath,jsonAnat);
        end
    end

end

