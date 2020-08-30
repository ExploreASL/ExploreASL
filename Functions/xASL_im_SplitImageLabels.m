function xASL_im_SplitImageLabels(ImagePaths, LabelTable, OutputFolder, bOverwrite)
%xASL_im_SplitImageLabels Extract individual label/regions from image(s)
%
% FORMAT: xASL_im_SplitImageLabels(ImagePaths, Path2TSV, OutputFolder)
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
%
% OUTPUT: n/a
% OUTPUT FILES: separate label images
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function allows extracting of labels from a NIfTI file
%              containing multiple labels, into single NIfTI files each
%              containing a single label.
%              Not all existing labels need to be extracted.
%              The following steps are performed:
%              1. Load TSV file
%              2. Process images
%
% EXAMPLE LabelTable as separate TSV-file: xASL_im_SplitImageLabels(xASL_adm_GetFileList(x.D.PopDir,'^4V_(?!MAP).*\.nii$', 'FPList'), '/ExampleTerritoryLabels.tsv', x.D.PopDir);
% EXAMPLE LabelTable as cell: xASL_im_SplitImageLabels(xASL_adm_GetFileList(x.D.PopDir,'^4V_(?!MAP).*\.nii$', 'FPList'), {1, 'ICA-L'; 2 'ICA-R'; 3, 'POS-L'; 4, 'POS-R'}, x.D.PopDir);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL
    

if nargin<4 || isempty(bOverwrite)
    bOverwrite = false;
end
if nargin<3 || isempty(OutputFolder)
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
        end
        
        if bResetOutputFolder
            OutputFolder = [];
        end
    end
end        
fprintf('\n');


end