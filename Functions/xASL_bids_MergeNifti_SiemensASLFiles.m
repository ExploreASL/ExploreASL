function pathOut = xASL_bids_MergeNifti_SiemensASLFiles(NiftiPaths)
%xASL_bids_MergeNifti_SiemensASLFiles Take a list of NIfTI files and
%concatenates 3D/4D files into a 4D sequence if possible (Siemens)
%
% FORMAT: pathOut = xASL_bids_MergeNifti_SiemensASLFiles(NiftiPaths)
% 
% INPUT:
%   NiftiPaths - cell containing list of strings with full paths of the files (REQUIRED)
%
% OUTPUT:
%   pathOut    - output paths
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Take a list of NIfTI files and concatenates 3D/4D files into a 4D sequence if possible (Siemens).
%
% EXAMPLE:     pathOut = xASL_bids_MergeNifti_SiemensASLFiles(NiftiPaths);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% xASL_bids_MergeNifti_SiemensASLFiles Merge Siemens ASL files with specific filename pattern

    bCheckConsistency = 1; % So far, no error was found and files can be concatenated
    bAlternatingControlLabel = 0; % == 1 when two files are given, one file is filled with controls and the other with labels
    pathOut = ''; % Newly assigned path of a concatenated file
    listEndNumber = zeros(length(NiftiPaths),1);

    for iFile=1:length(NiftiPaths)
        [Fpath, Ffile] = xASL_fileparts(NiftiPaths{iFile});
        jsonParms = spm_jsonread(fullfile(Fpath,[Ffile,'.json']));

        % Make sure that we deal with Siemens files
        if isempty(jsonParms) || ~isfield(jsonParms,'Manufacturer') || isempty(strfind(jsonParms.Manufacturer,'Siemens'))
            bCheckConsistency = 0;
        end

        % Get the size of the 4th dimension
        dimTime = size(xASL_io_Nifti2Im(NiftiPaths{iFile}),4);

        % Compare the PEAxis and check that it is the same for all of the M0-files
        if iFile == 1
            checkTimeDim = dimTime;
        else
            if checkTimeDim ~= dimTime
                bCheckConsistency = 0;
            end
        end
    end

    % Conditions passed for merging Siemens files
    if bCheckConsistency
        if checkTimeDim == 2
            % Alternating scenario C/L - merge them in the order of the last number in the filename
            for iFile=1:length(NiftiPaths)
                [Fpath, Ffile] = xASL_fileparts(NiftiPaths{iFile});

                % List the end number from the file name
                [iStart, iEnd] = regexp(Ffile,'\d*$');
                listEndNumber(iFile) = str2double(Ffile(iStart:iEnd));
            end
            [~,indexSortedFile] = sort(listEndNumber);

        elseif checkTimeDim == 1
            % Each volume contains only a label and control
            if length(NiftiPaths) == 2
                % If there are only two of them then sort by the trailing number in the filename
                for iFile=1:length(NiftiPaths)
                    [Fpath, Ffile, ~] = xASL_fileparts(NiftiPaths{iFile});

                    % List the end number from the file name
                    [iStart, iEnd] = regexp(Ffile,'\d*$');
                    listEndNumber(iFile) = str2double(Ffile(iStart:iEnd));
                end
                [~, indexSortedFile] = sort(listEndNumber);
            else
                % If there are multiple, then break to tags ASL4D_x_x_Y.nii and controls ASL4D_Y.nii
                % Sort each by Y and alternate them
                % If there are only two of them then sort by the trailing number in the filename
                listTag = zeros(length(NiftiPaths),1);
                for iFile=1:length(NiftiPaths)
                    [Fpath, Ffile, ~] = xASL_fileparts(NiftiPaths{iFile});

                    % List the end number from the file name
                    [iStart, iEnd] = regexp(Ffile,'\d*$');
                    listEndNumber(iFile) = str2double(Ffile(iStart:iEnd));

                    % Check for the x_x_Y pattern in the filename identifying the tags
                    [iStart, iEnd] = regexp(Ffile,'\d*_\d*_\d*$');
                    if ~isempty(iStart) && ~isempty(iEnd)
                        listTag(iFile) = 1;
                    else
                        listTag(iFile) = 0;
                    end

                end
                [~, indexSortedFileTag] = sort(listEndNumber);
                indexSortedFileTag = indexSortedFileTag(listTag(indexSortedFileTag) == 1);
                [~, indexSortedFileControl] = sort(listEndNumber);
                indexSortedFileControl = indexSortedFileControl(listTag(indexSortedFileControl) == 0);
                if length(indexSortedFileTag) == length(indexSortedFileControl)
                    indexSortedFile(1:2:length(NiftiPaths)) = indexSortedFileTag;
                    indexSortedFile(2:2:length(NiftiPaths)) = indexSortedFileControl;
                else
                    bCheckConsistency = 0;
                end
            end
        elseif checkTimeDim > 2 && length(NiftiPaths) == 2
            % Exactly two files with longer time dimension is given - we count that one is control and the other label and we merge them
            % by alternating the volumes from both
            indexSortedFile = [1 2];
            bAlternatingControlLabel = 1;
        else
            % The Siemens specific scenarios work only with single or two volumes per file
            % And also with the case that exactly two files are given with more than a single time-dim
            % All others are skipped and can be merged according to the generic criteria.
            bCheckConsistency = 0;
        end
    end

    % Correctly identified one of the Siemens scenarios so can merge the files in the indexSortedFile order
    if bCheckConsistency
        pathOut = xASL_bids_MergeNifti_Merge(NiftiPaths,indexSortedFile,'ASL4D',bAlternatingControlLabel);

        if ~isempty(pathOut)
            xASL_bids_MergeNifti_RenameParms(Fpath,'ASL4D');
            xASL_bids_MergeNifti_Delete(NiftiPaths);
            fprintf('Corrected dcm2niiX output for Siemens files:\n');
            fprintf('%s\n', pathOut);
        end
    end


end
