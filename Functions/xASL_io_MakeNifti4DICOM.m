function xASL_io_MakeNifti4DICOM(PathIn, x, DataType)
%xASL_io_MakeNifti4DICOM Create NIfTI that is ready for conversion to uint16 DICOM
%
% FORMAT: xASL_io_MakeNifti4DICOM(PathIn, x)
%
% INPUT:
%   PathIn - path to NIfTI file to convert (REQUIRED)
%   x      - struct containing pipeline environment parameters (REQUIRED)
%   Datatype = 'INT16' or 'UINT16' (OPTIONAL, DEFAULT = 'INT16')
%
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function converts a NIfTI file to a DICOM file for
%              PACS visualization purposes:
%              For scaling/visualization:
%              1) Remove peak signal
%              2) Remove valley signal
%              3) Remove NaNs
%              4) Rescale to 12 bit integers
%              5) Save NIfTI. We also zip the NIfTI as this NIfTI won't be opened by ExploreASL
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

CBFim = xASL_io_Nifti2Im(PathIn);
[Fpath, Ffile] = xASL_fileparts(PathIn);

FileNew = [Ffile '_Visual2DICOM.nii'];
RegExpNew = [Ffile '_Visual2DICOM\.nii'];
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

%% Save NIfTI
xASL_io_SaveNifti(PathIn, PathNew, CBFim, [], false);

newNifti = xASL_io_ReadNifti(PathNew);

if strcmp(DataType,'UINT16')
    newNifti.dat(:,:,:) = uint16(CBFim);
    newNifti.dat.scl_slope = 2/RescaleSlope;
    newNifti.dat.dtype = 'UINT16-LE';
else
    newNifti.dat(:,:,:) = int16(CBFim);
    newNifti.dat.scl_slope = 1/RescaleSlope;
    newNifti.dat.dtype = 'INT16-LE';
end
    
create(newNifti);
xASL_adm_ZipFileList(Fpath, RegExpNew, false, true, [0 Inf], true);


end