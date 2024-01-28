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
% EXAMPLE used by xASL_wrp_ResampleASL: xASL_io_ASLSubtractionAveraging(x, saveWhichNifti, 0, PathASL4D{iSpace});
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

% pathNifti     = {path_PWI               path_PWI3D               path_PWI4D                 path_Control                path_Control3D                  path_Control4D};

% [fPath, fFile]  = fileparts(x.P.Path_mean_control);
% x.P.Path_Control3D = fullfile(fPath, [fFile '3D.nii']);
% x.P.Path_Control4D = fullfile(fPath, [fFile '4D.nii']);
% 
% pathDefault   = {x.P.Path_PWI           x.P.Path_PWI3D           x.P.Path_PWI4D             x.P.Path_mean_control       x.P.Path_Control3D              x.P.Path_Control4D              x.P.Path_ASL4D};


%% ========================================================================================
%% 2. Load NIfTIs & JSONs
for iNifti=1:5
    if nargin<(iNifti+3)
        path2load{iNifti} = [];
    else
        path2load{iNifti} = varargin{iNifti};
    end
    [x, imageInput{iNifti}] = xASL_io_ASLSubtractionAveraging_sub_LoadPWI_JSON(x, path2load{iNifti});
end


%% ========================================================================================
%% 3. Run the subtraction/averaging (==ASL volume manipulation)
% imageInput = {ASL4D PWI4D PWI3D Control4D Control3D}
% imageOutput = {PWI PWI3D PWI4D Control Control3D Control4D}

% [PWI, PWI3D, PWI4D, x, Control, Control3D, Control4D] = xASL_im_ASLSubtractionAveraging(x, ASL4D [, PWI4D, PWI3D, Control4D, Control3D])
[imageOutput{1}, imageOutput{2}, imageOutput{3}, x, imageOutput{4}, imageOutput{5}, imageOutput{6}] = xASL_im_ASLSubtractionAveraging(x, imageInput{1}, imageInput{2}, imageInput{3}, imageInput{4}, imageInput{5});
 

%% 3. =====================================================================================
%% 4. Save the NIfTI & JSON files
fieldNamesPLD = {'InitialPLD_PWI'       'InitialPLD_PWI3D'       'InitialPLD_PWI4D'         'InitialPLD_Control'        'InitialPLD_Control3D'          'InitialPLD_Control4D'};
fieldNamesLD  = {'LabelingDuration_PWI' 'LabelingDuration_PWI3D' 'LabelingDuration_PWI4D'   'LabelingDuration_Control'  'LabelingDuration_Control3D'    'LabelingDuration_Control4D'};
fieldNamesTE  = {'EchoTime_PWI'         'EchoTime_PWI3D'         'EchoTime_PWI4D'           'EchoTime_Control'          'EchoTime_Control3D'            'EchoTime_Control4D'};

for iNifti=1:size(saveWhichNifti, 1)
    n2save = saveWhichNifti{iNifti,1};
    xASL_io_ASLSubtractionAveraging_sub_SavePWI_JSON(x, saveWhichNifti{iNifti,2}, imageOutput{n2save}, fieldNamesPLD{n2save}, fieldNamesLD{n2save}, fieldNamesTE{n2save}, bCopyOrigJson);
end


end




%% ========================================================================================
%% ========================================================================================
function [x, imageOut] = xASL_io_ASLSubtractionAveraging_sub_LoadPWI_JSON(x, path_Nifti)
%xASL_io_ASLSubtractionAveraging_sub_LoadJSON Load NIfTI image matrix & sidecar JSON

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
            jsonFields = fields(json.Q);
            fieldPLD = find(contains(jsonFields, 'Initial_PLD'));
            fieldLD = find(contains(jsonFields, 'LabelingDuration'));
            fieldTE = find(contains(jsonFields, 'EchoTime'));

            if ~isempty(fieldPLD) && ~isempty(fieldLD)
                x.Q.LabelingDuration = json.Q.LabelingDuration;
                x.Q.Initial_PLD = json.Q.Initial_PLD;

                if ~isempty(fieldTE)
                    x.Q.EchoTime = json.Q.EchoTime;
                else
                    x.Q.EchoTime = [];
                end
            end
        end
    end


end



%% ========================================================================================
%% ========================================================================================
function xASL_io_ASLSubtractionAveraging_sub_SavePWI_JSON(x, path2save, image2save, fieldNamePLD, fieldNameLD, fieldNameTE, bCopyOrigJson)
%xASL_io_ASLSubtractionAveraging_sub_LoadJSON Load NIfTI image matrix & sidecar JSON

    if ~isempty(image2save)
        % Also save parameter vectors
        jsonFields.Initial_PLD = x.Q.(fieldNamePLD);
        jsonFields.LabelingDuration = x.Q.(fieldNameLD);
    
        if isfield(x.Q, fieldNameTE)
            jsonFields.EchoTime = x.Q.(fieldNameTE);
        end
    
        xASL_io_SaveNifti(x.P.Path_ASL4D, path2save, single(image2save), 32, 0, [], bCopyOrigJson, xASL_bids_parms2BIDS(jsonFields, [], 1));
    end

end