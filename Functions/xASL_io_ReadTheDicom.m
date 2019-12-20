function [Info] = xASL_io_ReadTheDicom(bUseDCMTK, DicomPath)
%xASL_io_ReadTheDicom Wrapper around DICOMTK & dicominfo for reading the dicom fields
%
% FORMAT: [Info] = xASL_io_ReadTheDicom(bUseDCMTK, DicomPath)
% 
% INPUT:
%   bUseDCMTK   - true if we prefer to read DICOMs using DCMTK compilation, false if we prefer to use dicominfo. 
%                 DCMTK is faster whereas dicominfo is more flexible
%                 if false: we read the DICOM ith dicominfo
%                 if true: we try reading the DICOM with DCMTK, and if this
%                 fails (defined by not hving read the EchoTime, RepetitionTime
%                 or ImageType), we repeat reading the DICOM by DCMTK
% DicomPath     - path to the DICOM we aim to read
%
% OUTPUT: 
%   Info        - struct containing DICOM fields and their information
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function tries to read a DICOM and throws a warning if it fails to
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [Info] = xASL_io_ReadTheDicom(bUseDCMTK, DicomPath)
% __________________________________
% Copyright 2015-2019 ExploreASL


warning('off','images:dicominfo:fileVRDoesNotMatchDictionary');

try
    if bUseDCMTK
        Info = xASL_io_DcmtkRead(DicomPath);

        if isempty(Info.EchoTime) || isempty(Info.RepetitionTime) || isempty(Info.ImageType)
            Info = dicominfo(DicomPath);
            if ~isfield(Info,'EchoTime') || ~isfield(Info,'RepetitionTime') || ~isfield(Info,'ImageType')
                error('Both DCMTK output & dicominfo were incomplete');
            else
                error('DCMTK output was incomplete, but dicominfo was complete');
            end
        end
    else
        Info = dicominfo(DicomPath);
    end
catch ME
    warning(ME.message);
end

end