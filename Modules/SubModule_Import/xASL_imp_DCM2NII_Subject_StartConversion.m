function [imPar, globalCounts, x, summary_line, destdir, scanpath, scan_name, dcm2niiCatchedErrors, nii_files, first_match] = xASL_imp_DCM2NII_Subject_StartConversion(imPar, globalCounts, x, bSkipThisOne, summary_line, destdir, scanpath, scan_name, dcm2niiCatchedErrors)
%xASL_imp_DCM2NII_Subject_StartConversion Start of DCM2NII subject conversion.
%
% FORMAT: [imPar, globalCounts, x, summary_line, destdir, scanpath, scan_name, dcm2niiCatchedErrors, nii_files, first_match] = xASL_imp_DCM2NII_Subject_StartConversion(imPar, globalCounts, x, bSkipThisOne, summary_line, destdir, scanpath, scan_name, dcm2niiCatchedErrors)
%
% INPUT:
%   imPar                 - Structure with import parameters (REQUIRED, STRUCT)
%   globalCounts          - Converted, skipped & missing scans (REQUIRED, STRUCT)
%   x                     - ExploreASL x structure (REQUIRED, STRUCT)
%   nii_files             - NIfTI files
%   bSkipThisOne          - Skip this one (BOOLEAN, REQUIRED)
%   summary_line          - Summary line
%   destdir               - Destination directory
%   scanpath              - Scan path
%   scan_name             - Scan name
%   dcm2niiCatchedErrors  - DCM2NII catched errors
%
% OUTPUT:
%   Almost the same as input + nii_files & first_match
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Start of DCM2NII subject conversion.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     n/a
%
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Start the conversion if this scan should not be skipped
    if bSkipThisOne
        summary_line = sprintf(',"skipped",,,,,,,,');
        globalCounts.skipped_scans(iSubject, iVisit, iSession, iScan) = 1;
    else
        nii_files = {};
        xASL_adm_CreateDir(destdir);

        %% Check if we have a nii(gz) file, or something that needs to be converted (parrec/dicom)
        if ~exist(scanpath, 'dir') && ~isempty(regexpi(scanpath,'(\.nii|\.nii\.gz)$'))
            %% We found a NIfTI file, check if output exists
            first_match = fullfile(destdir, [scan_name '.nii']);
            if imPar.bOverwrite || ~xASL_exist(first_match,'file')
                [~, fname, fext] = fileparts(scanpath);
                destfile = fullfile(destdir, [fname fext]); % will be renamed later
                xASL_Copy(scanpath, destfile, imPar.bOverwrite, imPar.bVerbose);
                % gunzip if required
                destfile = xASL_adm_UnzipNifti(destfile);
                xASL_Move(destfile, first_match, imPar.bOverwrite, imPar.bVerbose);
            end
            nii_files{1} = first_match;
            
        else % we found dicom files
            %% Start the conversion. Note that the dicom filter is only in effect when a directory is specified as input.
            try
                [nii_files, scan_name, first_match, MsgDcm2nii] = xASL_io_dcm2nii(scanpath, destdir, scan_name, ...
                    'DicomFilter', imPar.dcmExtFilter, 'Verbose', imPar.bVerbose, 'Overwrite', imPar.bOverwrite, ...
                    'Version', imPar.dcm2nii_version, 'x', x);

                % If dcm2nii produced a warning or error, catch this & store it
                if ~isempty(MsgDcm2nii) && ~isempty(regexpi(MsgDcm2nii,'.*(error).*')) % if it contains a warning/error
                    dcm2niiCatchedErrors = xASL_imp_CatchErrors('xASL_io_dcm2nii', MsgDcm2nii, dbstack, ...
                        ['dcm2nii_' imPar.dcm2nii_version], pwd, scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar);
                end

            catch ME
                dcm2niiCatchedErrors = xASL_imp_CatchErrors(ME.identifier, ME.message, [], ...
                    [], [], scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar, ME.stack);

                if imPar.bVerbose; warning(['dcm2nii ' scanpath ' crashed, skipping']); end
                if imPar.bVerbose; warning('Check whether the scan is complete'); end
                first_match = xASL_adm_GetFileList(scanpath, ['.*' imPar.dcmExtFilter],'FPList',[0 Inf]);
                if  ~isempty(first_match); first_match = first_match{1}; end
            end
        end
    end

end


