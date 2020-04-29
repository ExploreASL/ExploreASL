function xASL_io_MakeNifti4DICOM(PathIn, x, DataType, OrientationPath, ResamplePath)
%xASL_io_MakeNifti4DICOM Create NIfTI that is ready for conversion to uint16 DICOM
%
% FORMAT: xASL_io_MakeNifti4DICOM(PathIn, x)
%
% INPUT:
%   PathIn - path to NIfTI file to convert (REQUIRED)
%   x      - struct containing pipeline environment parameters (REQUIRED)
%   Datatype = 'INT16' or 'UINT16' (OPTIONAL, DEFAULT = 'INT16')
%   OrientationPath = Path of NIfTI to align with (OPTIONAL, DEFAULT = empty, no registration)
%   ResamplePath = Path of NIfTI to resample to (OPTIONAL, DEFAULT = empty, no resampling)
%
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function converts a NIfTI file to one that is ready to convert to DICOM for
%              PACS visualization purposes:
%              For scaling/visualization:
%              1) Remove peak signal
%              2) Remove valley signal
%              3) Remove NaNs
%              4) Rescale to 12 bit integers
%              5) Save NIfTI. We also zip the NIfTI as this NIfTI won't be opened by ExploreASL
%              6) Manage scale slope/datatype
%              7) Apply original orientation
%              8) Zip NIfTI
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_io_MakeNifti4DICOM(x.P.PathCBF, x);
% __________________________________
% Copyright 2015-2019 ExploreASL


%% Admin
if nargin<3 || isempty(DataType)
    DataType = 'INT16'; % default
elseif isempty(regexp(DataType,'^(U|)INT16$'))
    error('Wrong DataType input parameter');
end
if nargin<4 || isempty(OrientationPath)
    OrientationPath = [];
    % by default put the CBF image back in the original T1w space
end
if nargin<5
    ResamplePath = []; % by default don't resample
end

CBFim = xASL_io_Nifti2Im(PathIn);
[Fpath, Ffile] = xASL_fileparts(PathIn);

FileNew = [Ffile '_Visual2DICOM'];
RegExpNew = [Ffile '_Visual2DICOM'];    

if ~isempty(OrientationPath)
    [~, FfileRef] = xASL_fileparts(OrientationPath);
    FileNew = [FileNew '_align-' FfileRef];
    RegExpNew = [RegExpNew '_align-' FfileRef];
end

if ~isempty(ResamplePath)
    [~, FfileRef] = xASL_fileparts(ResamplePath);
    FileNew = [FileNew '_space-' FfileRef];
    RegExpNew = [RegExpNew '_space-' FfileRef];
end

FileNew = [FileNew '.nii'];
RegExpNew = [RegExpNew '\.nii'];

PathNew = fullfile(Fpath, FileNew);


%% 1) Remove peak signal
[~, CBFim] = xASL_im_MaskPeakVascularSignal(CBFim, [], true, 6);

%% 2) Remove valley signal
% Define the lower threshold
NegativeSignal = CBFim(CBFim<0);

% medianValue = xASL_stat_MedianNan(NegativeSignal);
madValue = xASL_stat_MadNan(NegativeSignal,1);
ClipThr = 2*madValue; % medianValue - (4*madValue);

CBFim(CBFim<ClipThr) = 0;
CBFim(CBFim<0) = 0;

%% 3) Remove NaNs
CBFim(isnan(CBFim)) = 0;

%% 4) Rescale to 12 bit integers
RescaleSlope = 4095/max(CBFim(:));
CBFim = round(CBFim.*RescaleSlope);

%% 5) Save NIfTI
xASL_io_SaveNifti(PathIn, PathNew, CBFim, [], false);
newNifti = xASL_io_ReadNifti(PathNew);

%% 6) Manage scale slope/datatype
if strcmp(DataType,'UINT16')
    newNifti.dat(:,:,:) = uint16(CBFim);
    newNifti.dat.scl_slope = 1/RescaleSlope;
    newNifti.dat.dtype = 'UINT16-LE';
else
    newNifti.dat(:,:,:) = int16(CBFim);
    newNifti.dat.scl_slope = 1/RescaleSlope;
    newNifti.dat.dtype = 'INT16-LE';
end

%% 7) Apply original orientation
% Here, we assume that the ASL scan is registered with the T1 NIfTI.
% The original position of the T1 scan, can be found in the mat0 of the
% NIfTI, either the T1_ORI or the T1.

% PM: If the T1.nii doesnt exist, we can compare ASL4D.nii.mat with
% ASL4D.nii.mat0 instead, but this will revert any registrations of ASL4D
% with structural scans

% We assume that the pipeline has done the following:
% 1) original T1 orientation from scanner = T1[_ORI].nii.mat0
% 2) T1 (& derivatives) moved to MNI ACPC center
% 3) ASL (& derivatives) registered to T1
% -> Hence, we only revert step 1-2 here

if ~isempty(OrientationPath)
    % define path to backup NIfTI
    [Fpath, Ffile, Fext] = xASL_fileparts(OrientationPath);
    OrientationPathOri = fullfile(Fpath, [Ffile '_ORI.nii']);
    if ~xASL_exist(OrientationPathOri,'file')
        OrientationPathOri = OrientationPath;
    end
    
    
    niiReg = xASL_io_ReadNifti(OrientationPath);
    niiReg_ORI = xASL_io_ReadNifti(OrientationPathOri);
    matReg = niiReg.mat;
    matReg_ORI = niiReg_ORI.mat0;

    Transformation = matReg_ORI/matReg;

    % Apply this transformation (which should be translations & rotations only)
    newNifti.mat = Transformation*newNifti.mat;
end

create(newNifti);

%% 8) Resample
if ~isempty(ResamplePath)
    if ~xASL_exist(ResamplePath, 'file') || isempty(regexp(ResamplePath,'\.nii(|\.gz)$'))
        fprintf('Warning xASL_io_MakeNifti4DICOM: Tried to resample but this NIfTI file was missing or invalid:\n');
        fprintf('%s\n', ResamplePath);
    else
        % First create a temporary copy with the original orientation
        [Fpath3, Ffile3] = xASL_fileparts(ResamplePath);
        TempPath = fullfile(Fpath3, [Ffile3 '_temp.nii']);
        xASL_Copy(ResamplePath, TempPath, 1);
        nii = xASL_io_ReadNifti(TempPath);
        nii.mat = nii.mat0;
        create(nii);
        
        xASL_spm_reslice(TempPath, PathNew, [], [], x.Quality, PathNew);
        xASL_delete(TempPath);
    end
end        


%% 9) Zip NIfTI
xASL_adm_ZipFileList(Fpath, RegExpNew, false, true, [0 Inf], true);


end