function xASL_spm_deformations(x, PathIn, PathOut, Interpolation, InverseSpace, AffineTrans, DeformationPath, VoxelSize)
% xASL_spm_deformations ExploreASL wrapper for SPM deformation tool
% xASL_spm_deformations ExploreASL wrapper for SPM deformation tool
%
% FORMAT: xASL_spm_deformations([x,] PathIn, PathOut[, Interpolation, InverseSpace, AffineTrans, DeformationPath])
% 
% INPUT:
%   PathIn              - path or cell containing list of paths of images to transform/warp (REQUIRED)
%   PathOut             - path to cell containing list of paths of output names. Needs to be same numel as PathIn (REQUIRED)
%   x                   - ExploreASL structure containing fields with global information about the pipeline environment
%                         and settings (e.g. x.Quality), useful when you want this script to copy the options of an ExploreASL pipeline run (OPTIONAL)
%   Interpolation       - interpolation setting used by warping, options:
%                         0) Nearest Neighbor, 1) Trilinear, 2-7) 2nd-7th degree B-spline (OPTIONAL, DEFAULT=2)
%   InverseSpace        - path to space that you want to warp to, inversily applying the transformation field (e.g. when warping an MNI map to
%                         subject space) (OPTIONAL, DEFAULT=empty, no inverse transformation)
%   AffineTrans         - ASL affine transformation to add , e.g. if geometric distortion was corrected by a fieldmap or affine
%                         registration, here the distortion correction between ASL and T1w can be added (OPTIONAL, DEFAULT=empty)
%   DeformationPath     - if other deformation field needs to be applied (OPTIONAL, DEFAULT=empty)
%                         By default the '/Mypath/Subject1/y_T1.nii' transformation from the structural pipeline is used
%   VoxelSize           - if resulting files need to be resampled to a
%                         voxel size that differs from the ExploreASL
%                         default [1.5 1.5 1.5] mm (OPTIONAL, DEFAULT=[1.5 1.5 1.5])
%                         
% OUTPUT: n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This ExploreASL wrapper manages the SPM deformation tool.
% It takes multiple (ExploreASL pipeline) transformations and combines/concatenates them  
% into a single transformation prior to applying it to the input images. 
% This allows to apply multiple transformations with a single interpolation, avoiding
% propagation of undesired interpolation effects. Mainly used to get native
% space images into standard space, or vice versa.
% Best to combine as many files as possible within this function, since the
% deformation calculation (which is the most computation intensive part) needs to be performed once for multi-file resampling
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_spm_deformations(x, '/MyStudy/Subject1/T1.nii.gz', '/MyStudy/Population/rT1.nii.gz');
% __________________________________
% Copyright 2015-2019 ExploreASL

%% ------------------------------------------------------------------------------------------------------------
%% Administration
if (iscell(PathIn) && ~iscell(PathOut)) || (~iscell(PathIn) && iscell(PathOut))
    error('PathIn & PathOut should both be either cell/string & equally large!');
end

if nargin<1 || isempty(x)
    x.Quality = true; % default quality
elseif ~isfield(x,'Quality')
    x.Quality = true;
end
if ~isfield(x,'P')
    x.P = struct;
end
if ~isfield(x,'D')
    x.D = struct;
end
if ~isfield(x.P,'SubjectID')
    x.P.SubjectID = 'Unknown subject';
end
if ~isfield(x.P,'STRUCT')
    x.P.STRUCT = 'T1';
end
if ~isfield(x.D,'ROOT')
    x.D.ROOT = '';
end
if ~isfield(x.D,'PopDir')
    x.D.PopDir = fullfile(x.D.ROOT,'dartel');
end
if ~iscell(PathIn) % force cell format
    t_IN = PathIn;
    t_OUT = PathOut;

    clear PathIn PathOut
    PathIn{1} = t_IN;
    PathOut{1} = t_OUT;
end
if length(PathIn)~=length(PathOut)
    error('PathIn & PathOut should have the same length!');
end

if nargin<5 || isempty(InverseSpace)
    AddInverse = false;
elseif ~xASL_exist(InverseSpace,'file')
    warning([InverseSpace ' didnt exist, skipping...']);
    return;
else
    [~, InverseSpaceFile] = fileparts(InverseSpace);
    xASL_adm_UnzipNifti(InverseSpace); % inverse space is always within data folder, so unzipping "on site"
    AddInverse = true;
end
AddAffine = false;
if nargin<6 || isempty(AffineTrans)
    AffineTrans = '';
elseif ~xASL_exist(AffineTrans, 'file')
    [~, AffineFile, AffineExt] = fileparts(AffineTrans);
    fprintf([AffineFile AffineExt ' didnt exist, skipping...\n']);
else
    [~, AffineTransFile] = fileparts(AffineTrans);
    AddAffine = true;
end
if nargin<7 || isempty(DeformationPath)
    if ~strcmp(x.P.SubjectID, 'Unknown subject')
        DeformationPath = fullfile(x.D.ROOT, x.P.SubjectID, ['y_' x.P.STRUCT '.nii']); % defaults at this
    else
        warning('DeformationPath input parameter missing');
    end
end
if nargin<8 || isempty(VoxelSize)
    VoxelSize = [1.5 1.5 1.5];
end
    
if ~xASL_exist(DeformationPath, 'file')
    if isfield(x.P,'Path_y_ASL') && strcmp(DeformationPath, x.P.Path_y_ASL)
        xASL_wrp_CreateASLDeformationField(x); % backwards compatibility
    else
        warning([DeformationPath ' didnt exist, skipping...']);
        return;
    end
end

xASL_adm_UnzipNifti(DeformationPath); % always within data folder, can be unzipped "on site"





%% ------------------------------------------------------------------------------------------------------------
%% Interpolation admin
if nargin<4 || isempty(Interpolation)
    Interpolation = 2; % default interpolation
end
if ~x.Quality && Interpolation~=0
    Interpolation = 1; % faster, for lower quality (unless this was 0 before)
elseif Interpolation<0 || Interpolation>7
    error(['Wrong interpolation set: ' num2str(Interpolation)]);
end

InterpMeth = {'Nearest Neighbor' 'Trilinear' '2nd B-spline' '3rd degree Bspline' '4th degree Bspline' '5th degree Bspline' '6th degree Bspline' '7th degree Bspline'};

matlabbatch{1}.spm.util.defs.out{1}.pull.interp = Interpolation;

fprintf('\n%s\n','------------------------------------------------------------------------------------------');
fprintf('%s\n',['SPM deformation utility executed, using ' InterpMeth{Interpolation+1} ' interpolation']);


%% ------------------------------------------------------------------------------------------------------------
%% Skip unexisting files
tIN = '';
tOUT = '';
for iS=1:length(PathIn)

    if ~xASL_exist(PathIn{iS},'file')
        fprintf('%s\n',['Transforming ' PathIn{iS} ' skipped, was non-existing']);
    else
        fprintf('%s\n',['Transforming ' PathIn{iS} ' to ' PathOut{iS}]);
        tIN{end+1,1} = PathIn{iS};
        tOUT{end+1,1} = PathOut{iS};
    end
end

if isempty(tIN)
    fprintf('%s\n','Transforming skipped, no files entered');
    return;
end

PathIn = tIN;
PathOut = tOUT;

fprintf('%s',' ');

%% Remove .gz extensions
if strcmp(DeformationPath(end-2:end),'.gz')
    DeformationPath = DeformationPath(1:end-3);
end


%% ------------------------------------------------------------------------------------------------------------
%% If input & output directory are not the same, first copy the files to the output directory. & unzip
% & rename variables accordingly. This allows for using read-only
% files/directories, that can be zipped

for iS=1:length(PathIn)
    Fpath1 = fileparts(PathIn{iS});
    Fpath2 = fileparts(PathOut{iS});

    if ~strcmp(Fpath1,Fpath2) % if the directory names are not the same
        % 1) Delete PathOut if exists (otherwise this will go wrong with
        % copying/unzipping)
        xASL_delete(PathOut{iS});
        % 2) Copy INPUT file to OUTPUT file (this would be overwritten anyway, below)
        xASL_Copy(PathIn{iS},PathOut{iS});
        % 3) Rename the inputname variable accordingly
        PathIn{iS} = PathOut{iS};
    end
    xASL_adm_UnzipNifti(PathIn{iS}); % unzip if needed    
end


% Here, we determine what kind of joint transformation should be applied
fprintf('%s\n','Combining transformations:');
matlabbatch{1}.spm.util.defs.comp  = ''; % initialize transformation combination


%% If inverse transformation and x.Quality==0, then invert the transformation on lowest resolution
%  This means we will reslice the InverseSpaceFile to the y_Transformation
%  resolution, if the y_Transformation resolution is lower than the destination
%  resolution

RevertResolution = false;
if AddInverse && (Interpolation==0 || x.Quality==0)
    TransformationResolution = xASL_io_ReadNifti(DeformationPath);
    TransformationResolution = TransformationResolution.hdr.pixdim(2:4);
    InverseSpaceResolution = xASL_io_ReadNifti(InverseSpace);
    InverseSpaceResolution = InverseSpaceResolution.hdr.pixdim(2:4);
    if prod(TransformationResolution)>prod(InverseSpaceResolution)
        % * if the destination space has a higher resolution, resample the InverseSpace
        % * we do this with lowest quality, since we only need this NIfTI
        % orientation matrix as reference, we don't care for this NIfTI image
        % * we omit the affine transformation here (taken care of below)
        PathTempRefSpace = fullfile(fileparts(InverseSpace),'TempRefSpace.nii');
        xASL_spm_reslice(DeformationPath, InverseSpace, [], [], false, PathTempRefSpace, 0);
        InverseSpaceOriginal = InverseSpace;
        InverseSpace = PathTempRefSpace;
        RevertResolution = true;
    end
end
   
%% ------------------------------------------------------------------------------------------------------------
%% If an affine registration of ASL to T1w exists, we apply this
if AddAffine
    if AddInverse
        [~, AffineTransFile] = fileparts(AffineTrans);
        matlabbatch{1}.spm.util.defs.comp{end+1}.inv.comp{1}.sn2def.matname     = {AffineTrans};
        matlabbatch{1}.spm.util.defs.comp{end}.inv.comp{1}.sn2def.vox           = [NaN NaN NaN];
        matlabbatch{1}.spm.util.defs.comp{end}.inv.comp{1}.sn2def.bb            = [NaN NaN NaN
                                                                                   NaN NaN NaN]; 
        matlabbatch{1}.spm.util.defs.comp{end}.inv.space                        = {x.D.ResliceRef};
        fprintf('%s\n',['Inverse affine transformation ' AffineTransFile '.mat, to space: ' InverseSpaceFile '.nii']);
    else
        matlabbatch{1}.spm.util.defs.comp{end+1}.sn2def.matname     = {AffineTrans};
        matlabbatch{1}.spm.util.defs.comp{end}.sn2def.vox           = [NaN NaN NaN];
        matlabbatch{1}.spm.util.defs.comp{end}.sn2def.bb            = [NaN NaN NaN
                                                                       NaN NaN NaN];
        fprintf('%s\n',['Affine transformation ' AffineTransFile '.mat']);                                                                       
    end
end


%% ------------------------------------------------------------------------------------------------------------
%% Here, we apply the transformation to MNI (inversely if requested)

if AddInverse
    matlabbatch{1}.spm.util.defs.comp{end+1}.inv.comp{1}.def    = {DeformationPath};
    matlabbatch{1}.spm.util.defs.comp{end}.inv.space            = {InverseSpace};
    fprintf('%s\n',['Inverse segmentation transformation y_' x.P.STRUCT '.nii of ' x.P.SubjectID ', to space: ' InverseSpaceFile '.nii']);
else
    fprintf('%s\n',['Segmentation transformation y_' x.P.STRUCT '.nii of ' x.P.SubjectID]);
    matlabbatch{1}.spm.util.defs.comp{end+1}.def                = {DeformationPath};
end

%% ------------------------------------------------------------------------------------------------------------
%% Determine the output voxelsize
if min(VoxelSize~=[1.5 1.5 1.5])
    matlabbatch{1}.spm.util.defs.comp{end+1}.idbbvox.vox = VoxelSize;
    matlabbatch{1}.spm.util.defs.comp{end}.idbbvox.bb = [NaN NaN NaN
                                                       NaN NaN NaN];
end

%% ------------------------------------------------------------------------------------------------------------
%% Define path to apply deformations to
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = PathIn;
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1; % saves in SUBJECTDIR, because we need to rename still
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];

%% ------------------------------------------------------------------------------------------------------------
%% Run the job & rename
spm_jobman('run',matlabbatch); % this applies the SPM job (i.e. joint transformation)

for iL=1:length(PathIn)
    [Fpath, Ffile, Fext] = fileparts(PathIn{iL});
    wname = fullfile(Fpath, ['w' Ffile Fext]);
    if xASL_exist(PathIn{iL}, 'file') && ~strcmp(wname, PathOut{iL})
        xASL_Move(wname, PathOut{iL}, 1);
    end
end

if RevertResolution
    for iS=1:length(PathOut)
        xASL_spm_reslice(InverseSpaceOriginal, PathOut{iS}, [], [], x.Quality, PathOut{iS}, Interpolation);
    end
end

if exist('PathTempRefSpace','var') % if we created a temporary space to calculate inverse transformation to
    xASL_delete(PathTempRefSpace);
end


end
