function xASL_imp_NII2BIDS_Subject(imPar, bidsPar, studyPar, nameSubject)
%xASL_imp_NII2BIDS_Subject Run NII to ASL-BIDS for one individual subject.
%
% FORMAT: xASL_imp_NII2BIDS_Subject(imPar, bidsPar, studyPar, nameSubject)
% 
% INPUT:
%   imPar           - JSON file with structure with import parameter (REQUIRED, STRUCT)
%   bidsPar         - Output of xASL_imp_Config (REQUIRED, STRUCT)
%   studyPar        - JSON file with the BIDS parameters relevant for the whole study (REQUIRED, STRUCT)
%   nameSubject     - name of the subject (REQUIRED, LIST)
%
% OUTPUT:
%   n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Run NII to ASL-BIDS for one individual subject.
%
% 1. Initialize
% 2. Process all the anatomical files
% 3. Process the perfusion files (iterate over sessions)
% - 1. Make a subject directory
% - 2. Iterate over runs
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     xASL_imp_NII2BIDS_Subject(imPar, bidsPar, studyPar, nameSubject);
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% 1. Initialize
    
    subjectLabel = xASL_adm_CorrectName(nameSubject,2);

	% Make a subject directory
	if ~exist(fullfile(imPar.BidsRoot,['sub-' subjectLabel]),'dir')
		mkdir(fullfile(imPar.BidsRoot,['sub-' subjectLabel]));
	end

    %% 2. Process all the anatomical files
    % Go throught the list of anat files
    for iAnatType = bidsPar.listAnatTypes

        % Check if it exists
        anatPath = '';
        if xASL_exist(fullfile(imPar.TempRoot,nameSubject,[iAnatType{1},'.nii']),'file')
            anatPath = fullfile(imPar.TempRoot,nameSubject,iAnatType{1});
        end

        if xASL_exist(fullfile(imPar.TempRoot,nameSubject,[iAnatType{1} '_1'],[iAnatType{1},'.nii']),'file')
            anatPath = fullfile(imPar.TempRoot,nameSubject,[iAnatType{1} '_1'],iAnatType{1});
        end

        % If anatomical file of this type exist, then BIDSify its structure
        if ~isempty(anatPath)

            % Create the anatomical directory
            if ~exist(fullfile(imPar.BidsRoot,['sub-' subjectLabel],'anat'),'dir')
                mkdir(fullfile(imPar.BidsRoot,['sub-' subjectLabel],'anat'));
            end

            % Move the NiFTI file
            xASL_Move([anatPath '.nii'],fullfile(imPar.BidsRoot,['sub-' subjectLabel],'anat',...
                ['sub-' subjectLabel '_' iAnatType{1} '.nii.gz']),1);

            % Load the JSON
            jsonAnat = spm_jsonread([anatPath,'.json']);
			
            % Save the JSON
			jsonAnat = xASL_bids_BIDSifyAnatJSON(jsonAnat,studyPar);
            jsonAnat = xASL_bids_VendorFieldCheck(jsonAnat);
            jsonAnat = xASL_bids_JsonCheck(jsonAnat,'');
			
            spm_jsonwrite(fullfile(imPar.BidsRoot,['sub-' subjectLabel],'anat',['sub-' subjectLabel '_' iAnatType{1} '.json']),jsonAnat);
        end
    end

    %% 3. Process the perfusion files
    listSessions = xASL_adm_GetFileList(fullfile(imPar.TempRoot,nameSubject),'^ASL.+$',false,[],true);

    % Go through all sessions
    for iSession = 1:length(listSessions)
        xASL_imp_NII2BIDS_SubjectSession(imPar, bidsPar, studyPar, listSessions, nameSubject, subjectLabel, iSession);
    end
end
