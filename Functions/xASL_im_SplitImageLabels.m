function xASL_im_SplitImageLabels(ImagePaths, LabelTable, OutputFolder, bOverwrite, ResampleDir, SubRegExp)
%xASL_im_SplitImageLabels Extract individual label/regions from image(s)
%
% FORMAT: xASL_im_SplitImageLabels(ImagePaths, LabelTable[, OutputFolder, bOverwrite, ResampleDir, SubRegExp])
%
% INPUT:
%   ImagePaths      - cell containing paths to label images (REQUIRED)
%   LabelTable      - cell containing label numbers and names for 
%                     each label that needs to be extracted.
%                     Can also be a path to a CSV or TSV file containing this cell.
%                     (REQUIRED)
%   OutputFolder    - path to folder where all separate label images will be stored
%                     when empty, each file will be outputted in the folder
%                     of the original image (OPTIONAL, DEFAULT = empty)
%   bOverwrite      - true for overwriting existing NIfTIs
%                     (OPTIONAL, DEFAULT = false)
%   ResampleDir     - path to folder where resampled label image to standard space
%                     should be saved. If empty, resampling is skipped. (OPTIONAL, DEFAULT = empty)
%   SubRegExp       - regular expression of subject name/ID inside the
%                     filepath, used when resampling to standard space
%                     e.g. x.subject_regexp (OPTIONAL, DEFAULT = empty).
%
% OUTPUT: 
%  variables        - n/a
%  files            - separate label images
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% DESCRIPTION: This function allows extracting of labels from a NIfTI file
%              containing multiple labels, into single NIfTI files each
%              containing a single label.
%              Not all existing labels need to be extracted.
%
%              The following steps are performed:
%
%              1. Load TSV file
%              2. Process images
%
% EXAMPLE: LabelTable as separate TSV-file: xASL_im_SplitImageLabels(xASL_adm_GetFileList(x.D.PopDir,'^4V_(?!MAP).*\.nii$', 'FPList'), '/ExampleTerritoryLabels.tsv', x.D.PopDir);
%          LabelTable as cell: xASL_im_SplitImageLabels(xASL_adm_GetFileList(x.D.ROOT,'^4V\.nii$', 'FPListRec'), {1, 'ICA-L'; 2 'ICA-R'; 3, 'POS-L'; 4, 'POS-R'}, [], 1, x.D.PopDir, x.subject_regexp);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% Copyright 2015-2020 ExploreASL
    

if nargin<6 || isnumeric(SubRegExp)
    SubRegExp = [];
end
if nargin<5 || isnumeric(ResampleDir)
    ResampleDir = [];
elseif isempty(SubRegExp)
    error('Resampling requested but no subject regular expression provided');
end
if nargin<4 || isempty(bOverwrite)
    bOverwrite = false;
end
if nargin<3 || isempty(OutputFolder) || isnumeric(OutputFolder)
    warning('No output folder specified, outputting file in its original folder');
    OutputFolder = [];
end
if nargin<2 || isempty(LabelTable)
    error('Please provide the label table as cell or path to a CSV/TSV file');
end
if nargin<1 || isempty(ImagePaths)
    error('Please provide (a list of) path(s) to image(s) containing label(s)');
end

if ~iscell(ImagePaths)
    ImagePaths = {ImagePaths};
end

if isempty(OutputFolder)
    bResetOutputFolder = true;
else
    bResetOutputFolder = false;
end

%% -------------------------------
%% 1. Load label table
if ~iscell(LabelTable) && ischar(LabelTable)
    [~, ~, Fext] = xASL_fileparts(LabelTable);

    if strcmp(Fext, '.csv')
        try
            LabelTable = xASL_csvRead(LabelTable);
        catch ME
            warning('CSV-file detected, but couldnt load it');
            error(ME);
        end        
    elseif strcmp(Fext, '.tsv')
        try
            LabelTable = xASL_tsvRead(LabelTable);
        catch ME
            warning('TSV-file detected, but couldnt load it');
            error(ME);
        end
    else
        error('Unrecognized file format, should have CSV or TSV extension');
    end
end
    
if ~iscell(LabelTable)
    error('Unrecognized label table format, should be a cell or path');
else
    try
        fprintf('Loading label table\n');
        for iLabel=1:size(LabelTable,1)
            LabelNumber(iLabel) = xASL_str2num(LabelTable{iLabel,1});
            LabelName{iLabel} = LabelTable{iLabel,2};
        end
    catch ME
        warning('Couldnt load label table, something wrong with the format?');
        error(ME);
    end
end
    

%% -------------------------------
%% 2. Process images
fprintf('Process images:   ');

for iImage=1:length(ImagePaths)
    xASL_TrackProgress(iImage, length(ImagePaths));
    % Get original filename
    [Fpath, Ffile] = xASL_fileparts(ImagePaths{iImage});
    
    % Load image
    ThisImage = xASL_io_Nifti2Im(ImagePaths{iImage});
    
    % Loop through labels
    for iLabel = 1:length(LabelNumber)
        LabelImage = ThisImage==LabelNumber(iLabel);
        if isempty(OutputFolder)
            OutputFolder = Fpath;
        end
        FileName = fullfile(OutputFolder, [LabelName{iLabel} '_' Ffile '.nii']);
        if xASL_exist(FileName, 'file') && ~bOverwrite
            fprintf('%s\n', [FileName ' already existed, skipping']);
        else
            xASL_io_SaveNifti(ImagePaths{iImage}, FileName, LabelImage);
            
            if ~isempty(ResampleDir)
                % Get subject id
                [StartIndex, EndIndex] = regexp(FileName, SubRegExp(2:end-1));
                SubjectID = FileName(StartIndex:EndIndex);
                % get filetype
                [~, FileType] = xASL_fileparts(FileName);
                
                PathOut = fullfile(ResampleDir, [FileType '_' SubjectID '.nii']);
                Path_y_ASL = fullfile(FileName(1:EndIndex), 'ASL_1', 'y_ASL.nii');
                if ~xASL_exist(Path_y_ASL)
                    warning(['file didnt exist: ' Path_y_ASL]);
                elseif ~xASL_exist(FileName)
                    warning(['file didnt exist: ' FileName]);
                else
                    xASL_spm_deformations([], FileName, [PathOut '.gz'], 0, [], [], Path_y_ASL);
                end
            end
        end
        
        if bResetOutputFolder
            OutputFolder = [];
        end
    end
end        
fprintf('\n');


end