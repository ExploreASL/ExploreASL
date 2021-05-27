function [x] = xASL_stat_AtlasForStats(x)
%xASL_stat_AtlasForStats Load atlases for ROI statistics
%
% FORMAT: [x] = xASL_stat_AtlasForStats(x)
%
% INPUT:
%   x                   - struct containing statistical pipeline environment parameters (REQUIRED)
%   x.S.InputAtlasPath  - path to NIfTI file containing atlas to load (REQUIRED)
%   x.utils.WBmask      - WholeBrain mask used to convert image to column &
%                         vice versa (ExploreASL compression method) (REQUIRED)
%   x.S.ROInamesPath    - path to TSV file containing ROI names for atlas.
%                         This TSV file should contain 1 row, with each
%                         cell corresponding to the ROI number (i.e. an
%                         atlas with 6 ROIs should have a TSV file with 6
%                         cells/columns) (OPTIONAL, DEFAULT=ask for input
%                         OR use ROI_1 ROI_2 ROI_n as names)
%   x.S.SubjectWiseVisualization - true to evaluate ROI masks for each subject,
%                                  takes computational time (OPTIONAL, DEFAULT=false)
%
% OUTPUT:
%   x                   - same as input
%   x.S.InputMasks      - ROI masks to compute statistics for, converted/compressed to columns
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function loads atlases, checks them, and
%              their ROI names, for later use as ROI definition in xASL_stat_GetROIstatistics
%              Note that the atlases should be integer values, or different 4rd
%              dimensions (i.e. multiple images), that are mutually
%              exclusive. This function takes the following steps:
%
%              1. Load atlas ROI names
%                 There should be a TSV sidecar to the atlas NIfTI file, as
%                 explained above.
%              2. deal with memory mapping
%              3. Resample atlas 50 1.5 mm^3 MNI
%              4. Converted atlas with integers to 4D binary image
%              5. Convert/compress masks into Columns
%              6. Print atlas overview image
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: x = xASL_stat_AtlasForStats(x);
% __________________________________
% Copyright 2015-2019 ExploreASL


%% Admin
% Check different file types
[Fpath, Ffile, Fext]  = xASL_fileparts(x.S.InputAtlasPath);
if ~strcmp(Fext,'dat') && isempty(findstr(Fext,'nii'))
    if exist(fullfile(Fpath,[Ffile '.dat']),'file')
           x.S.InputAtlasPath = fullfile(Fpath,[Ffile '.dat']);
    elseif xASL_exist(fullfile(Fpath,[Ffile '.nii']))
         x.S.InputAtlasPath = fullfile(Fpath,[Ffile '.nii']);
    else
         warning('incorrect atlas selected, please try again');
         fprintf('%s\n', x.S.InputAtlasPath);
         return;
    end
elseif ~xASL_exist(x.S.InputAtlasPath, 'file')
     warning('incorrect atlas selected, please try again');
     fprintf('%s\n', x.S.InputAtlasPath);
     return;
end

SumMask = sum(x.utils.WBmask(:));

if strcmp(x.S.InputAtlasPath(end-2:end),'.gz')
    x.S.InputAtlasPath = x.S.InputAtlasPath(1:end-3); % allow .gz or .nii input
end

%% 1) Load atlas ROI names
x.S.ROInamesPath = fullfile(Fpath,[Ffile '.tsv']);

% Find out whether we are run in GUI or CLI
 if ~usejava('desktop') || ~usejava('jvm') || ~feature('ShowFigureWindows')
     UseGUI = false;
 else
     UseGUI = true;
 end

if ~exist(x.S.ROInamesPath,'file')
    if UseGUI
        YesNo = questdlg(['Would you like to specify a tsv-file containing the ROI names for ' x.S.InputAtlasPath '?']);
    else
        fprintf('%s\n',['Would you like to specify a tsv-file containing the ROI names for ' x.S.InputAtlasPath '?']);
        fprintf('%s\n','Otherwise we will use ROI_1 ROI_2 ROI_n as names of the atlas ROI masks');
        YesNo = input('Please type Yes or No');
    end
    if strcmp(YesNo,'Yes')
        if UseGUI
            [Fname, Fpathstr] = uigetfile('*.tsv', ['Select the tab-separated value file containing the ROI names for ' x.S.InputAtlasPath]);
            if sum(Fpathstr==0) || ~strcmp(Fname(end-3:end),'.tsv')
                return
            else
                x.S.ROInamesPath = fullfile(Fpathstr,Fname);
            end
        else
            x.S.ROInamesPath = input('Please enter the path to the TSV file containing the ROI names:');
        end
    end
end

if isfield(x.S,'ROInamesPath') && exist(x.S.ROInamesPath, 'file') % open TSV file
    fclose all;
    FileID = fopen(x.S.ROInamesPath);

    C = textscan(FileID,'%s');
    C = C{1};

    if ~isempty(C)
        % Concatenate cells, if multiple cells
        if length(C)>1
            Cnew = C{1};
            for iC=2:length(C)
                Cnew =[Cnew ' ' C{iC}];
            end
        else
            Cnew = C{1};
        end

        CommasIndex = strfind(Cnew,' ');

        if isempty(CommasIndex)
            x.S.NamesROI{1} = Cnew;
        else
            x.S.NamesROI{1} = Cnew(1:CommasIndex(1)-1);
            for iC=2:length(CommasIndex) % work around for missing strsplit in previous matlab versions
                x.S.NamesROI{iC} = Cnew(CommasIndex(iC-1)+1:CommasIndex(iC)-1);
            end
            if ~strcmp(Cnew(end),',')
                iC = length(CommasIndex);
                x.S.NamesROI{iC+1} = Cnew(CommasIndex(iC)+1:end);
            end
        end
    end
    fclose(FileID);
end


%% 2) Load atlas image matrix, deal with memory mapping
AtlasIsColumns = false;
if ischar(x.S.InputAtlasPath) % allows both image input or ImagePath input
    if strcmp(x.S.InputAtlasPath(end-3:end),'.dat') % if memory mapping, then load this
        %% Part for Atlas stored as columns
        AtlasIsColumns = true;
        TempAtlas = memmapfile(x.S.InputAtlasPath);

        if isfield(x.S,'NamesROI') % get number of masks
            nMasks = length(x.S.NamesROI);
        else
            nMasks = 6;
        end
        nSubj = size(TempAtlas.Data,1)/SumMask/nMasks;
        InputAtlasIM = reshape(TempAtlas.Data,[SumMask nMasks nSubj]); % reshape into [Brainvoxels nMasks nSubjects]
    else
        InputAtlasIM = xASL_io_Nifti2Im(x.S.InputAtlasPath);
    end
else
    InputAtlasIM = x.S.InputAtlasPath;
end

if ~AtlasIsColumns
    %% 3) Resample atlas 50 1.5 mm^3 MNI
    % Allow for multiple voxel dimensions, BUT THIS FORCES MNI 1.5 mm
    SizeAtlas = size(InputAtlasIM);
    DimRatioAtlas = SizeAtlas(1:3)./[121 145 121];
    if prod(DimRatioAtlas)~=1
        fprintf('%s\n','Atlas has different dimensions, make sure it is in MNI space');
        if round(DimRatioAtlas(1),2)==round(DimRatioAtlas(2),2) && round(DimRatioAtlas(2),2)==round(DimRatioAtlas(3),2)
            fprintf('%s\n','Detected that atlas has different MNI dimensions, resampling to 1.5 mm');
            [Fpath, Ffile, Fext] = xASL_fileparts(x.S.InputAtlasPath);
            Atlas15Path = fullfile(Fpath,[Ffile '1.5mm' Fext]);
            xASL_spm_reslice( x.D.ResliceRef, x.S.InputAtlasPath, [], [], x.Quality, Atlas15Path, 0 );
        end
        InputAtlasIM = xASL_io_Nifti2Im(Atlas15Path);
    end

    %% 4) Converted atlas with integers to 4D binary image
    %  Allow for multiple atlas forms (3D or 4D), later transformed to multi-atlas 2D (Columns)
    if ~(size(InputAtlasIM,4)==1 && max(InputAtlasIM(:))>1)
        % don't need to reorganize
    else
        AtlasIn = InputAtlasIM;
        AtlasOut = zeros([size(AtlasIn(:,:,:,1,1,1)) max(AtlasIn(:))],'uint8');
        for iL=1:max(AtlasIn(:))
            tempIM = zeros(size(AtlasIn(:,:,:,1,1,1)));
            tempIM(AtlasIn==iL) = 1;
            AtlasOut(:,:,:,iL) = tempIM;
        end
        InputAtlasIM = AtlasOut;
    end

    %% 5) Convert/compress masks into Columns
    fprintf('%s\n','Converting masks:   ');
    x.S.InputMasks = zeros(sum(x.utils.WBmask(:)),size(InputAtlasIM,4),'uint8'); % memory pre-allocation
    for iL=1:size(InputAtlasIM,4)
        xASL_TrackProgress(iL,size(InputAtlasIM,4));
        x.S.InputMasks(:,iL,:) = xASL_im_IM2Column(InputAtlasIM(:,:,:,iL,[1:size(InputAtlasIM,5)]),x.utils.WBmask);
    end
    fprintf('\n');
else
    x.S.InputMasks = InputAtlasIM;
end

%% Create dummy ROI names, if we don't have them
if ~isfield(x.S,'NamesROI')
    if size(InputAtlasIM,4)>1
        maxROI = size(InputAtlasIM,4);
    else
        maxROI = max(InputAtlasIM(:));
    end
    for iR=1:maxROI
        x.S.NamesROI{iR} = ['ROI_' num2str(iR)]; % default ROIs
    end
end

%% 6) Print atlas overview image (takes time, disabled by default)
if x.S.SubjectWiseVisualization
    fprintf('Printing subject-specific masks (if exist) together in label colors:   ')
    % CAVE: only one label per voxel will be shown (latest have preference,
    % can be improved later)
    for iSub=1:size(x.S.InputMasks,3)
        xASL_TrackProgress(iSub,size(x.S.InputMasks,3));
        LabelIM = xASL_vis_TransformData2View(xASL_Convert4D_3D_atlas(xASL_im_Column2IM(x.S.InputMasks(:,:,iSub),x.utils.WBmask)),x);
        DataIM = xASL_vis_TransformData2View(x.utils.skull.*xASL_io_Nifti2Im(fullfile(x.D.TemplateDir,'rT1.nii')),x);
        CombiIM = xASL_im_ProjectLabelsOverData(DataIM,LabelIM,x);
        xASL_vis_Imwrite(CombiIM, fullfile(x.S.CheckMasksDir,[Ffile '_Subj_' num2str(iSub) '.jpg']));
    end
end


end
