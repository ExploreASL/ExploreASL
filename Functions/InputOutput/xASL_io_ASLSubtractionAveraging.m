function xASL_io_ASLSubtractionAveraging(x, saveWhichNifti, bCopyOrigJson, varargin)
%xASL_io_ASLSubtractionAveraging Wrapper for xASL_io_ASLSubtractionAveraging, managing the loading and saving
%
% FORMAT: xASL_io_ASLSubtractionAveraging(x, saveWhichNifti, bCopyOrigJson, path_ASL4D, path_PWI4D, path_PWI3D, path_Control4D, path_Control3D)
%
% INPUT:
%   x               - structure containing fields with all information required to run this function (REQUIRED)
%   saveWhichNifti  - cell array containing column 1) indices of NIfTIs to save & column 2) paths to save NIfTI to (REQUIRED)
%                     1 = PWI
%                     2 = PWI3D
%                     3 = PWI4D
%                     4 = Control
%                     5 = Control3D
%                     6 = Control4D
%                     e.g.: {1, 'path/PWI.nii'; 5, 'path/Control3D.nii'}
%                     this concerns the output
%   bCopyOrigJson   - boolean stating if the original JSON contents are copied as well
%                     (OPTIONAL, DEFAULT = true)
%                     the following varargin parameters concern the input:
%   path_ASL4D      - Path to ASL4D NIfTI (STRING, OPTIONAL)
%   path_PWI4D      - likewise, but then for PWI4D
%   path_PWI3D      - likewise, but then for PWI3D
%   path_Control4D  - likewise, but then for Control4D
%   path_Control3D  - likewise, but then for Control3D
%
% OUTPUT:             n/a (all output is written to NIfTI and JSON files)
%                     see saveWhichNifti for the NIfTI output options
%
% DESCRIPTION: This function wraps around xASL_io_ASLSubtractionAveraging, for loading and saving (see "help xASL_io_ASLSubtractionAveraging"),
%              performing the following steps:
% 1. Admin
% 2. Load NIfTIs & JSONs
% 3. Run the subtraction/averaging (==ASL volume manipulation)
% 4. Save the NIfTI & JSON files
%
% Nomenclature:
% ASL4D = original timeseries
% PWI4D = subtracted timeseries -> for quantification purposes (also applies to control4D)
% PWI3D = averaged per TE—PLD—LabDur combination -> for QC purposes  (also applies to control3D)
% PWI   = single perfusion-weighted volume -> for registration purposes  (also applies to control)
% 
% % we use the same nomenclature for control: control4D, control3D, control, as these are the same control-label pairs but without subtraction
%
% EXAMPLE: 
%     Inside xASL_wrp_ResampleASL: xASL_io_ASLSubtractionAveraging(x, saveWhichNifti, 0, PathASL4D{iSpace});
%     Inside xASL_wrp_ResampleASL: xASL_io_ASLSubtractionAveraging(x, {4, x.P.Path_mean_control}, 1, x.P.Path_rdespiked_ASL4D)
%     Inside xASL_wrp_RegisterASL: xASL_io_ASLSubtractionAveraging(x, {1, x.P.Path_mean_PWI_Clipped;4, x.P.Path_mean_control}, 0, x.P.Path_despiked_ASL4D);
% __________________________________
% Copyright (C) 2015-2024 ExploreASL

%% ========================================================================================
%% 1. Admin & checks
% Input image volumes should be correct
% Creation of the PWI3D or PWI4D image volumes can be skipped if they are provided as input

if nargin<1 || isempty(x)
    error('x input missing');
end
if nargin<2 || isempty(saveWhichNifti)
    error('Second input argument missing');
end
if nargin<3 || isempty(bCopyOrigJson)
    bCopyOrigJson = true;
end


%% ========================================================================================
%% 2. Load NIfTIs & JSONs
% ASL4D, PWI4D, PWI3D, Control4D, Control3D

fieldNamesPLD = {'Initial_PLD'      'InitialPLD_PWI4D'       'InitialPLD_PWI3D'       'InitialPLD_Control4D'       'InitialPLD_Control3D'};
fieldNamesLD  = {'LabelingDuration' 'LabelingDuration_PWI4D' 'LabelingDuration_PWI3D' 'LabelingDuration_Control4D' 'LabelingDuration_Control3D'};
fieldNamesTE  = {'EchoTime'         'EchoTime_PWI4D'         'EchoTime_PWI3D'         'EchoTime_Control4D'         'EchoTime_Control3D'};

for iNifti=1:5
    if nargin<(iNifti+3) % if the path is not defined, we will skip loading
        path2load{iNifti} = []; % empty path
    else
        path2load{iNifti} = varargin{iNifti}; % path to load
    end
    [x, imageInput{iNifti}] = xASL_io_ASLSubtractionAveraging_sub_LoadPWI_JSON(x, path2load{iNifti}, fieldNamesPLD{iNifti}, fieldNamesLD{iNifti}, fieldNamesTE{iNifti});
end


%% ========================================================================================
%% 3. Run the subtraction/averaging (==ASL volume manipulation)
% imageInput = {ASL4D PWI4D PWI3D Control4D Control3D}
% imageOutput = {PWI PWI3D PWI4D Control Control3D Control4D}

% [PWI, PWI3D, PWI4D, x, Control, Control3D, Control4D] = xASL_im_ASLSubtractionAveraging(x, ASL4D [, PWI4D, PWI3D, Control4D, Control3D]);
% While the images 2-5 can be passed is images, the first has to be passed as a path as ASL4D might need to have motion correction applied
[imageOutput{1}, imageOutput{2}, imageOutput{3}, x, imageOutput{4}, imageOutput{5}, imageOutput{6}] = xASL_im_ASLSubtractionAveraging(x, path2load{1}, imageInput{2}, imageInput{3}, imageInput{4}, imageInput{5});
 

%% 3. =====================================================================================
%% 4. Save the NIfTI & JSON files
fieldNamesPLD = {'InitialPLD_PWI'       'InitialPLD_PWI3D'       'InitialPLD_PWI4D'         'InitialPLD_Control'        'InitialPLD_Control3D'          'InitialPLD_Control4D'};
fieldNamesLD  = {'LabelingDuration_PWI' 'LabelingDuration_PWI3D' 'LabelingDuration_PWI4D'   'LabelingDuration_Control'  'LabelingDuration_Control3D'    'LabelingDuration_Control4D'};
fieldNamesTE  = {'EchoTime_PWI'         'EchoTime_PWI3D'         'EchoTime_PWI4D'           'EchoTime_Control'          'EchoTime_Control3D'            'EchoTime_Control4D'};

for iNifti=1:size(saveWhichNifti, 1)
    n2save = saveWhichNifti{iNifti,1}; % which NIfTI to save

	% For saving, we have to use the path of the closest reference image
	switch n2save
		case 1
			pathReferenceList = [1,2,3];	% Saving 1 = PWI, check path 1,2,3
		case 2
			pathReferenceList = [1,2]; % Saving 2 = PWI3D, check path 1,2
		case 3
			pathReferenceList = 1; % Saving 3 = PWI4D, check path 1
		case 4
			pathReferenceList = [1,4,5]; % Saving 4 = Control, check path 1,4,5
		case 5
			pathReferenceList = [1,4]; % Saving 5 = Control3D, check path 1,4
		case 6
			pathReferenceList = 1; % Saving 6 = Control4D, check path 1
	end

	% By default use Path_ASL4D as the reference
	pathReference = x.P.Path_ASL4D;
	for iReference = 1:length(pathReferenceList)
		if ~isempty(path2load{pathReferenceList(iReference)})
			pathReference = path2load{pathReferenceList(iReference)};
		end
	end
    xASL_io_ASLSubtractionAveraging_sub_SavePWI_JSON(x, saveWhichNifti{iNifti,2}, imageOutput{n2save}, fieldNamesPLD{n2save}, fieldNamesLD{n2save}, fieldNamesTE{n2save}, bCopyOrigJson, pathReference);
end


end




%% ========================================================================================
%% ========================================================================================
function [x, imageOut] = xASL_io_ASLSubtractionAveraging_sub_LoadPWI_JSON(x, path_Nifti, fieldNamePLD, fieldNameLD, fieldNameTE)
%xASL_io_ASLSubtractionAveraging_sub_LoadJSON Load NIfTI image matrix & sidecar JSON
%
%  x               - structure containing fields with all information required to run this function (REQUIRED)
% path_Nifti      - path to NIfTI file (STRING, REQUIRED)
% fieldNamePLD    - field name for PostLabelingDelay (NUMERIC, REQUIRED)
% fieldNameLD     - field name for LabelingDuration (NUMERIC, REQUIRED)
% fieldNameTE     - field name for EchoTime (NUMERIC, REQUIRED)
% imageOut        - NIfTI image matrix

    imageOut = []; % default

    %% Avoiding loading if NIfTI path is not defined
    if isempty(path_Nifti)
        return;
    elseif iscell(path_Nifti)
        path_Nifti = path_Nifti{1};
    end
    if ~xASL_exist(path_Nifti, 'file')
        warning(['Missing file: ' path_Nifti]);
        return;
    end

    %% Load NIfTI
    imageOut = xASL_io_Nifti2Im(path_Nifti);

    %% Load JSON
    [fPath, fFile] = xASL_fileparts(path_Nifti);
    pathJson = fullfile(fPath, [fFile '.json']);
	if ~exist(pathJson, 'file')
        warning(['Reverting to x.Q memory: JSON sidecar missing, : ' pathJson]);
    else
        json = xASL_io_ReadJson(pathJson); % load the json-file
        json = xASL_bids_parms2BIDS([], json, 0); % Convert it from BIDS to Legacy

        % if it has the two required quantification fields PLD & LD, we will use this
        if isfield(json, 'Q')
			if isfield(json.Q, 'LabelingDuration')
				x.Q.(fieldNameLD) = json.Q.LabelingDuration;
			else
				x.Q.(fieldNameLD) = [];
				if ~strcmp(x.Q.LabelingType, 'pasl') % Labeling duration is required for PCASL, but not for PASL, there are saturation times there instead
					warning(['LD missing in: '  pathJson]);
				end
			end

			if isfield(json.Q, 'Initial_PLD')
				x.Q.(fieldNamePLD) = json.Q.Initial_PLD;
			else
				warning(['PLD missing in: '  pathJson]);
			end

			if isfield(json.Q, 'EchoTime')
				x.Q.(fieldNameTE) = json.Q.EchoTime;
			else
				warning(['EchoTime missing from: '  pathJson]);
			end
        end
	end

end

%% ========================================================================================
%% ========================================================================================
function xASL_io_ASLSubtractionAveraging_sub_SavePWI_JSON(x, path2save, image2save, fieldNamePLD, fieldNameLD, fieldNameTE, bCopyOrigJson, pathReference)
%xASL_io_ASLSubtractionAveraging_sub_LoadJSON Load NIfTI image matrix & sidecar JSON
%
% x               - structure containing fields with all information required to run this function (REQUIRED)
% path2save       - path to save NIfTI file (STRING, REQUIRED)
% image2save      - NIfTI image matrix (REQUIRED)
% fieldNamePLD    - field name for PostLabeling Delay (NUMERIC, REQUIRED)
% fieldNameLD     - field name for LabelingDuration (NUMERIC, REQUIRED)
% fieldNameTE     - field name for EchoTime (NUMERIC, REQUIRED)
% bCopyOrigJson   - boolean stating if the original JSON contents are copied as well

    if ~isempty(image2save)
        % Also save parameter vectors
        jsonFields.Initial_PLD = x.Q.(fieldNamePLD);
        jsonFields.LabelingDuration = x.Q.(fieldNameLD);
    
        if isfield(x.Q, fieldNameTE)
            jsonFields.EchoTime = x.Q.(fieldNameTE);
        end
    
        xASL_io_SaveNifti(pathReference, path2save, image2save, 32, 0, [], bCopyOrigJson, jsonFields, true);
    end

end