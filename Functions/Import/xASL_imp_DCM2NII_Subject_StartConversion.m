function [globalCounts, x, summary_line, destdir, scanpath, scan_name, dcm2niiCatchedErrors, nii_files, first_match] = xASL_imp_DCM2NII_Subject_StartConversion(globalCounts, x, bSkipThisOne, summary_line, destdir, scanpath, scan_name, dcm2niiCatchedErrors, scanFields)
%xASL_imp_DCM2NII_Subject_StartConversion Start of DCM2NII subject conversion.
%
% FORMAT: [globalCounts, x, summary_line, destdir, scanpath, scan_name, dcm2niiCatchedErrors, nii_files, first_match] = xASL_imp_DCM2NII_Subject_StartConversion(globalCounts, x, bSkipThisOne, summary_line, destdir, scanpath, scan_name, dcm2niiCatchedErrors)
%
% INPUT:
%   globalCounts           - Converted, skipped & missing scans (REQUIRED, STRUCT)
%   x                      - ExploreASL x structure (REQUIRED, STRUCT)
%   x.modules.import.imPar - Structure with import parameters (REQUIRED, STRUCT)
%   bSkipThisOne           - Skip this one (BOOLEAN, REQUIRED)
%   summary_line           - Summary line (REQUIRED)
%   destdir                - Destination directory (CHAR ARRAY, PATH, REQUIRED)
%   scanpath               - Scan path (CHAR ARRAY, PATH, REQUIRED)
%   scan_name              - Scan name (CHAR ARRAY, REQUIRED)
%   dcm2niiCatchedErrors   - DCM2NII catched errors (STRUCT, REQUIRED)
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
% Copyright 2015-2023 ExploreASL

    %% Start the conversion if this scan should not be skipped
    first_match = [];
    if bSkipThisOne
        summary_line = sprintf(',"skipped",,,,,,,,');
        globalCounts.skipped_scans(scanFields.iSubject, scanFields.iVisit, scanFields.iSession, scanFields.iScan) = 1;
    else
        nii_files = {};
        xASL_adm_CreateDir(destdir);

        %% Check if we have a nii(gz) file, or something that needs to be converted (parrec/dicom)
        if ~exist(scanpath, 'dir') && ~isempty(regexpi(scanpath,'(\.nii|\.nii\.gz)$'))
            %% We found a NIfTI file, check if output exists
            first_match = fullfile(destdir, [scan_name '.nii.gz']);
            % only copy this file if it doesn't exist yet, or if we want to
            % overwrite
            if x.modules.import.imPar.bOverwrite || ~xASL_exist(first_match,'file') % this will check both .nii & .nii.gz
                xASL_Copy(scanpath, first_match, x.modules.import.imPar.bOverwrite, x.modules.import.imPar.bVerbose);
				
				[jsonDir, jsonFile, ~] = xASL_fileparts(scanpath);
				if xASL_exist(fullfile(jsonDir, [jsonFile '.json']), 'file')
					xASL_Copy(fullfile(jsonDir, [jsonFile '.json']), fullfile(destdir, [scan_name '.json']), x.modules.import.imPar.bOverwrite, x.modules.import.imPar.bVerbose);
				end
                % gzip if required
                xASL_adm_GzipAllFiles(fileparts(first_match));
            end
            nii_files{1} = first_match;
            
        else %% Convert DICOM files
            %% Start the conversion. Note that the dicom filter is only in effect when a directory is specified as input.
            try
                [nii_files, scan_name, first_match, MsgDcm2nii] = xASL_io_dcm2nii(scanpath, destdir, scan_name, x.modules.import.imPar, x.opts.MyPath);

                % If dcm2nii produced a warning or error, catch this & store it
                if ~isempty(MsgDcm2nii) && ~isempty(regexpi(MsgDcm2nii,'.*(error).*')) % if it contains a warning/error
                    dcm2niiCatchedErrors = xASL_imp_CatchErrors('xASL_io_dcm2nii', MsgDcm2nii, dbstack, ...
                        ['dcm2nii_' x.modules.import.imPar.dcm2nii_version], pwd, scan_name, scanpath, destdir, dcm2niiCatchedErrors, x.modules.import.imPar);
                end

            catch ME
                dcm2niiCatchedErrors = xASL_imp_CatchErrors(ME.identifier, ME.message, [], ...
                    [], [], scan_name, scanpath, destdir, dcm2niiCatchedErrors, x.modules.import.imPar, ME.stack);

                % Print warnings in verbose mode
                if x.modules.import.imPar.bVerbose
                    warning(['dcm2nii ' scanpath ' crashed, skipping...']);
                end
                if x.modules.import.imPar.bVerbose
                    warning('Check whether the scan is complete...');
                end
                first_match = xASL_adm_GetFileList(scanpath, ['.*' x.modules.import.imPar.dcmExtFilter],'FPList',[0 Inf]);
                if  ~isempty(first_match); first_match = first_match{1}; end
            end
        end
    end

end


