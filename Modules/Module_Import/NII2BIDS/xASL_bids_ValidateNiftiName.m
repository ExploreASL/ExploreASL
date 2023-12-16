function xASL_bids_ValidateNiftiName(fileName,perfType)
%xASL_bids_ValidateNiftiName Validate the NIFTI filename based on the regular expressions from bids-matlab
%
% FORMAT: xASL_bids_ValidateNiftiName(fileName,perfType)
%
% INPUT:
%   fileName   - Insert the output NIFTI path of the asl or m0scan (STRING REQUIRED)
%   perfType   - Define whether we have an asl or m0scan (STRING, REQUIRED)
%
% OUTPUT: 
%   n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:     Validate the NIFTI filename based on the regular expressions from bids-matlab.
%
% EXAMPLE:         n/a
%
% __________________________________
% Copyright 2015-2021 ExploreASL


    if strcmp(perfType,'asl')
        % ASL
        labels = regexp({fileName}, [ ...
            '^sub-[a-zA-Z0-9]+' ...              % sub-<participant_label>
            '(?<ses>_ses-[a-zA-Z0-9]+)?' ...     % ses-<label>
            '(?<acq>_acq-[a-zA-Z0-9]+)?' ...     % acq-<label>
            '(?<rec>_rec-[a-zA-Z0-9]+)?' ...     % rec-<label>
            '(?<dir>_dir-[a-zA-Z0-9]+)?' ...     % dir-<label>
            '(?<run>_run-[a-zA-Z0-9]+)?' ...     % run-<index>
            '_asl\.nii(\.gz)?$'], 'names'); % NIfTI file suffix/extension
    elseif strcmp(perfType,'m0scan')
        % M0
        labels = regexp({fileName}, [ ...
            '^sub-[a-zA-Z0-9]+' ...              % sub-<participant_label>
            '(?<ses>_ses-[a-zA-Z0-9]+)?' ...     % ses-<label>
            '(?<acq>_acq-[a-zA-Z0-9]+)?' ...     % acq-<label>
            '(?<rec>_rec-[a-zA-Z0-9]+)?' ...     % rec-<label>
            '(?<dir>_dir-[a-zA-Z0-9]+)?' ...     % dir-<label>
            '(?<run>_run-[a-zA-Z0-9]+)?' ...     % run-<index>
            '_m0scan\.nii(\.gz)?$'], 'names'); % NIfTI file suffix/extension
    else
        warning('Unknown perfusion file type...');
        labels = cell(1,1);
    end
    
    % Check labels
    if isempty(labels{1,1})
        warning('Invalid NIFTI filename %s...',fileName);
    end


end



