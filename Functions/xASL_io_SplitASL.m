function xASL_io_SplitASL(inPath, iM0, iDummy)
%xASL_io_SplitASL Splits ASL & M0 & Dummy images, when they are within the same NIfTI
%
% FORMAT: xASL_io_SplitASL(inPath[, iM0, iDummy])
%
% INPUT:
%   inPath      - path to ASL NIfTI file (e.g. //analysis/ASL_1/ASL4D.nii) (REQUIRED)
%   iM0         - index/indices of volume(s) containing the M0 (OPTIONAL, DEFAULT = [])
%   iDummy      - index/indices of volume(s) containing the Dummy scans (OPTIONAL, DEFAULT = [])
%
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function splits ASL4D & M0 & Dummy images if they were in the same sequence.
%              If dcm2niiX has already splitted the ASL4D NIfTI, this is reconstructed first.
%              If no M0 exists, or only ASL splitting is desired, leave iM0 empty ([]).
%              The dummy scans can be excluded from the ASL sequence during the splitting. Both iM0 and iDummy
%              are the absolute positions of both in the original time series
%
%              Vendor product sequence examples:
%              GE 3D spiral sometimes puts the M0 at the last volume of the series -> iM0 = [2];
%              Philips 3D GRASE puts the M0 as control-label volume pair -> iM0 = [1 2];
%              Siemens 3D GRASE puts the M0 as the first volume -> iM0 = 1;
%              Some Siemens 3D GRASE puts a second Dummy control image -> iDummy = 2;
%
% 1. Input parameter admin 
% 2. Prepare paths
% 3. First concatenate NIfTIs
% 4. Save M0 NIfTI
% 5. Determine ASL indices
% 6. Save ASL4D NIfTI
% 7. Split relevant JSON parameters/arrays
% 8. Split ASL4Dcontext.tsv
% 9. Modify JSON fields
% 10. Copy sidecars
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: for moving the first two volumes to M0.nii:
% xASL_io_SplitASL('/data/RAD/share/EPAD500/010EPAD00001/ASL_1/ASL4D.nii', [1 2]);
% EXAMPLE 2: for concatenating files only:
% xASL_io_SplitASL('/data/RAD/share/EPAD500/010EPAD00001/ASL_1/ASL4D.nii', []);
% EXAMPLE 3: for moving the first two volumes to M0.ii and removing the 3rd volume as dummy:
% xASL_io_SplitASL('/data/RAD/share/EPAD500/010EPAD00001/ASL_1/ASL4D.nii', [1 2],3);
% EXAMPLE 4: for moving the third volume to M0.ii and removing the first 2 volumes as dummy:
% xASL_io_SplitASL('/data/RAD/share/EPAD500/010EPAD00001/ASL_1/ASL4D.nii', 3,[1 2]);
% EXAMPLE 5: for removing the 1st volume as dummy
% xASL_io_SplitASL('/data/RAD/share/EPAD500/010EPAD00001/ASL_1/ASL4D.nii', [],1);
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% -----------------------------------------------------------------------------------------------------------------------------------------------------
    %% 1. Input parameter admin
	% Can be char or number on input, output should be a horizontal vector

	if nargin < 2 || isempty(iM0)
		iM0 = [];
	else
		iM0 = xASL_str2num(iM0);
		if size(iM0,1)>1 && size(iM0,2)>1
			warning('M0PositionInASL4D should be provided as a vector and not as a matrix')
		end
		iM0 = iM0(:)';
	end
	
	if nargin < 3 || isempty(iDummy)
		iDummy = [];
	else
		iDummy = xASL_str2num(iDummy);
		if size(iDummy,1)>1 && size(iDummy,2)>1
			warning('DummyScanPositionInASL4D should be provided as a vector and not as a matrix')
		end
		iDummy = iDummy(:)';
	end
	
    % Only run the splitting if indexes are provided
    if isempty(iM0) && isempty(iDummy)
        return;
    end
    
	%% -----------------------------------------------------------------------------------------------------------------------------------------------------
	%% 2. Prepare paths
    [Fpath, Ffile] = xASL_fileparts(inPath);
    [paths] = xASL_io_SplitASL_GetPaths(Fpath, Ffile);
    [ASLlist] = xASL_io_SplitASL_RestoringFromBackup(Fpath, Ffile, paths, iM0);
    
    if ~xASL_exist(paths.ASL_Source,'file') % otherwise was already split

        %% 3. First concatenate NIfTIs
        xASL_io_SplitASL_ConcatNiftis(paths, ASLlist);
        
		% Load the backup-file
		tIM = xASL_io_Nifti2Im(paths.ASL_Source);
		
        %% 4. Save M0 NIfTI
        if ~xASL_exist(paths.M0,'file') && ~isempty(iM0) % don't overwrite, and skip if no M0
            xASL_io_SaveNifti(paths.ASL_Source,paths.M0 ,tIM(:,:,:,iM0),[],0);
            bCreateM0 = true;
        else
            bCreateM0 = false;
        end

        %% 5. Determine ASL indices
        ASLindices = 1:1:size(tIM,4);
        IndexASL   = ones(1,size(tIM,4));
		if ~isempty(iM0)
			IndexASL(iM0) = 0;
		end
		if ~isempty(iDummy)
			IndexASL(iDummy) = 0;
		end
        ASLindices          = ASLindices(logical(IndexASL));
		
        %% 6. Save ASL4D NIfTI
        xASL_io_SaveNifti(paths.ASL_Source,paths.ASL,tIM(:,:,:,ASLindices),[],false);
        
        %% 7. Split relevant JSON parameters/arrays
        [jsonM0, jsonASL] = xASL_io_SplitASL_SplitJSON(paths.ASL_Source_JSON, iM0, ASLindices, iDummy);

        %% 8. Split ASL4Dcontext.tsv
		xASL_io_SplitASL_SplitASLContext(paths.ASL_Source_TSV, paths.ASLTSV, iM0, ASLindices, iDummy);
        
        %% 9. Modify JSON fields
        [jsonM0, jsonASL] = xASL_io_SplitASL_PostModify(paths.ASL, jsonM0, jsonASL, iM0);
        
        %% 10. Copy sidecars
        if exist(paths.ASL_Source_MAT,'file') && bCreateM0
            xASL_Copy(paths.ASL_Source_MAT, paths.M0MAT);
        end
        if exist(paths.ASL_Source_MAT,'file') && ~strcmp(paths.ASL_Source_MAT, paths.ASLMAT)
            xASL_Copy(paths.ASL_Source_MAT, paths.ASLMAT);
        end
        if exist('jsonM0','var') && bCreateM0
            % BIDS validation
            jsonM0 = xASL_bids_VendorFieldCheck(jsonM0);
            jsonM0 = xASL_bids_JsonCheck(jsonM0,'M0');
            % Write file
            spm_jsonwrite(paths.M0JSON, jsonM0);
        end
        if exist('jsonASL','var') && ~strcmp(paths.ASL_Source_JSON, paths.ASLJSON)
            spm_jsonwrite(paths.ASLJSON, jsonASL);
        end
        
    end

end


%% xASL_io_SplitASL_ConcatNiftis
function xASL_io_SplitASL_ConcatNiftis(paths, ASLlist)

    if length(ASLlist)>1
        % Reconstruct the ASL4D first
        for iNii=1:length(ASLlist)
            tnii = xASL_io_ReadNifti(ASLlist{iNii});
            RescaleSlope(iNii) = tnii.dat.scl_slope;

            tIM = xASL_io_Nifti2Im(tnii.dat(:,:,:,:,:,:,:));
            Dim4 = size(tIM,4);
            if iNii==1
                TotIm = tIM;
            else
                TotIm(:,:,:,CDim4+1:CDim4+Dim4) = tIM;
            end
            CDim4 = size(TotIm,4);
        end

        if length(unique(RescaleSlope))>1
            warning('Rescaleslopes were not the same between concatenated NIfTIs, skipping...');
            return;
        else
            % Save concatenated ASL series as backup ASL
            xASL_io_SaveNifti(ASLlist{1}, paths.ASL_Source, TotIm, [], false);
            %
            if exist(paths.originalMAT,'file') && ~strcmp(paths.originalMAT, paths.ASL_Source_MAT)
                xASL_Move(paths.originalMAT, paths.ASL_Source_MAT);
            end
            if exist(paths.originalJSON,'file') && ~strcmp(paths.originalJSON, paths.ASL_Source_JSON)
                xASL_Move(paths.originalJSON, paths.ASL_Source_JSON);
            end

            for iNii=1:length(ASLlist)
                xASL_delete(ASLlist{iNii});
                [Fpath, Ffile] = xASL_fileparts(ASLlist{iNii});
                JSON2Delete = fullfile(Fpath, [Ffile '.json']);
                MAT2Delete = fullfile(Fpath, [Ffile '_parms.mat']);
                xASL_delete(JSON2Delete);
                xASL_delete(MAT2Delete);
            end
        end

    else % backup the ASL4D.nii & sidecars
        xASL_Move(ASLlist{1}, paths.ASL_Source, true);
        if exist(paths.originalMAT,'file') && ~strcmp(paths.originalMAT, paths.ASL_Source_MAT)
            xASL_Move(paths.originalMAT, paths.ASL_Source_MAT);
        end
        if exist(paths.originalJSON,'file') && ~strcmp(paths.originalJSON, paths.ASL_Source_JSON)
            xASL_Move(paths.originalJSON, paths.ASL_Source_JSON);
        end
    end

end


%% xASL_io_SplitASL_RestoringFromBackup
function [ASLlist] = xASL_io_SplitASL_RestoringFromBackup(Fpath, Ffile, paths, iM0)

    ASLlist = xASL_adm_GetFileList(Fpath, ['^' Ffile '(|_.)'  '(|_\d)' '\.nii$'], 'FPList', [0 Inf]);
  
    % Restore NIFTIs & JSONs from backup
    if xASL_exist(paths.ASL_Source,'file')
        fprintf('Restoring ASL4D.nii from backup NIfTI ASL4D_Source.nii\n');
        ASLlist{1} = paths.ASL;
        xASL_Move(paths.ASL_Source, ASLlist{1}, true);
        % Restore json as well
        if exist(paths.ASLJSON,'file') && exist(paths.ASL_Source_JSON,'file')
            xASL_Move(paths.ASL_Source_JSON, paths.ASLJSON, true);
        end
        % Restore TSV from backup
        if xASL_exist(paths.ASL_Source_TSV,'file')
            xASL_Move(paths.ASL_Source_TSV, paths.ASLTSV, true);
        end
        % Delete M0, allowing to restore this from the backup in case M0 is in the ASL timeseries
        if xASL_exist(paths.M0,'file') && ~isempty(iM0)
            fprintf('Deleting M0 to restore it from ASL4D\n');
            xASL_delete(paths.M0);
            xASL_delete(paths.M0JSON);
        end
    elseif isempty(ASLlist)
        error([paths.ASL ' didnt exist, skipping']);
    end

end


%% xASL_io_SplitASL_GetPaths
function [paths] = xASL_io_SplitASL_GetPaths(Fpath, Ffile)

    % Path ASL & M0 NIFTI
    paths.M0 = fullfile(Fpath,'M0.nii');
    paths.ASL = fullfile(Fpath, 'ASL4D.nii');
    
    % JSON & MAT
    paths.originalJSON = fullfile(Fpath,[Ffile '.json']);
    paths.originalMAT = fullfile(Fpath,[Ffile '_parms.mat']);
    
    % ASL MAT, JSON & TSV
    paths.ASLMAT = fullfile(Fpath, 'ASL4D_parms.mat');
    paths.ASLJSON = fullfile(Fpath, 'ASL4D.json');
    paths.ASLTSV = fullfile(Fpath, 'ASL4Dcontext.tsv');
    
    % M0 MAT & JSON
    paths.M0MAT = fullfile(Fpath, 'M0_parms.mat');
    paths.M0JSON = fullfile(Fpath, 'M0.json');
    
    % Source files
    paths.ASL_Source = fullfile(Fpath,'ASL4D_Source.nii.gz');
    paths.ASL_Source_MAT = fullfile(Fpath, 'ASL4D_Source_parms.mat');
    paths.ASL_Source_JSON = fullfile(Fpath, 'ASL4D_Source.json');
    paths.ASL_Source_TSV = fullfile(Fpath, 'ASL4Dcontext_Source.tsv');

end


%% xASL_io_SplitASL_PostModify
function [jsonM0, jsonASL] = xASL_io_SplitASL_PostModify(ASLname, jsonM0, jsonASL, iM0)

    % Check that both structs are not empty dummy structs
    if ~isempty(fieldnames(jsonASL)) && ~isempty(fieldnames(jsonM0))
        
        % Add IntendedFor field
        [~, file, extension] = xASL_fileparts(ASLname);
        jsonM0.IntendedFor = ['perf/' file extension];
        
        % Remove M0Type field of M0 JSON if it still exists
        if isfield(jsonM0, 'M0Type')
            jsonM0 = rmfield(jsonM0, 'M0Type');
        end

        % Change M0Type field of ASL JSON to Separate
        if isfield(jsonASL, 'M0Type') && ~isempty(iM0)
            jsonASL.M0Type = 'Separate';
        end

    end

end


%% xASL_io_SplitASL_SplitJSON
function [jsonM0, jsonASL] = xASL_io_SplitASL_SplitJSON(BackupJSONPath, indicesM0, indicesASL, indicesDummy)

    % Make sure that we have column arrays
    if size(indicesM0,2)>size(indicesM0,1)
        indicesM0 = indicesM0';
    end
    if size(indicesASL,2)>size(indicesASL,1)
        indicesASL = indicesASL';
    end
    if size(indicesDummy,2)>size(indicesDummy,1)
        indicesDummy = indicesDummy';
    end

    % Load backup JSON
    if exist(BackupJSONPath,'file')
        jsonStruct = spm_jsonread(BackupJSONPath);
        
		% Remove DummyScanPositionInASL4D and M0PositionInASL4D fields if present
		if isfield(jsonStruct,'M0PositionInASL4D')
			jsonStruct = rmfield(jsonStruct,'M0PositionInASL4D');
		end
		if isfield(jsonStruct,'DummyScanPositionInASL4D')
			jsonStruct = rmfield(jsonStruct,'DummyScanPositionInASL4D');
		end
		
        % Define M0 & ASL fallback JSONs
        jsonM0 = jsonStruct;
        jsonASL = jsonStruct;
        
        % Fields which should be split
        splitFields = {'EchoTime', 'FlipAngle', 'RepetitionTimePreparation', 'PostLabelingDelay', 'VascularCrushingVENC', 'LabelingDuration'}';
        
        % Iterate over JSON fields
        jsonFields = fieldnames(jsonStruct);
        for iField = 1:numel(jsonFields)
            % Check if fieldname is in splitFields list
            if ismember(jsonFields{iField,1},splitFields)
                % Check if field has to be split
                if numel(jsonStruct.(jsonFields{iField,1}))>1
                    % Compare array length
                    if (numel(indicesM0)+numel(indicesASL)+numel(indicesDummy))==numel(jsonStruct.(jsonFields{iField,1}))
                        % Get current array
                        currentArray = jsonStruct.(jsonFields{iField,1});
                        % Add correct field value to M0 json
                        jsonM0.(jsonFields{iField,1}) = currentArray(indicesM0);
                        % Add correct field value to ASL json
                        jsonASL.(jsonFields{iField,1}) = currentArray(indicesASL);
                    end
                end
            end
        end
        
    else
        % Fallback
        jsonM0 = struct;
        jsonASL = struct;
    end

end

%% xASL_io_SplitASL_SplitJSON
function xASL_io_SplitASL_SplitASLContext(BackupTSVPath, newTSVPath, indicesM0, indicesASL, indicesDummy)

% Make sure that we have column arrays
if size(indicesM0,2)>size(indicesM0,1)
	indicesM0 = indicesM0';
end
if size(indicesASL,2)>size(indicesASL,1)
	indicesASL = indicesASL';
end
if size(indicesDummy,2)>size(indicesDummy,1)
	indicesDummy = indicesDummy';
end

% Backup the TSV file
if xASL_exist(newTSVPath,'file') && (~isempty(indicesM0) || ~isempty(indicesDummy))
	xASL_Move(newTSVPath, BackupTSVPath);
end

% Load backup JSON
if exist(BackupTSVPath,'file')
	% Read the original context file
	aslContext = xASL_tsvRead(BackupTSVPath);
        
	% Remove non-ASL fields, but keep the header
	aslContext = aslContext([1;indicesASL+1]);
	
	% Save the new ASL context
	xASL_tsvWrite(aslContext, newTSVPath, true);
	
	% Check that all controls and labels come in pairs
	aslContext = aslContext(2:end);

	bidsPar = xASL_bids_Config;
	nControls = numel(regexpi(strjoin(aslContext),bidsPar.stringControl));
	nLabels = numel(regexpi(strjoin(aslContext),bidsPar.stringLabel));
	if nControls ~= nLabels
		warning('After removing M0 volume and Dummy scans, non-matching number of volumes (control-label)');
	end
	
end

end

