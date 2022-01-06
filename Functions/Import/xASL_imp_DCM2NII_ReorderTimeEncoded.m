function xASL_imp_DCM2NII_ReorderTimeEncoded(nii_files, bTimeEncoded, resultJSON)
%xASL_imp_DCM2NII_ReorderTimeEncoded Reorder TEs and PLDs accordingly for time encoded sequences
%
% FORMAT: xASL_imp_DCM2NII_ReorderTimeEncoded(nii_files, bTimeEncoded, resultJSON)
% 
% INPUT:
%  nii_files     - List of nifti files (REQUIRED)
%  bTimeEncoded  - Time encoded sequence or not (BOOLEAN, REQUIRED)
%  resultJSON    - JSON structure (STURCT, REQUIRED)
%
% OUTPUT:
%   n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Reorder TEs and PLDs accordingly for time encoded sequences.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     n/a
%
% __________________________________
% Copyright 2015-2022 ExploreASL

    if numel(nii_files)>=1
        [resultPath, resultFile] = xASL_fileparts(nii_files{1});
        % Check if we the current sequence is a Hadamard or not
        if bTimeEncoded
            % Check image
            if xASL_exist(nii_files{1},'file')
                % Determine the number of time points within each NIfTI
                imASL = xASL_io_Nifti2Im(nii_files{1});
                numberTEs = length(resultJSON.EchoTime);
                numberPLDs = int32(size(imASL,4)/numberTEs);
                
                % Reorder TEs and PLDs - first cycle TE afterwards PLD
                vectorOldOrder = zeros(size(imASL,4),1);
                for iPLD = 1:(double(numberPLDs))
                    vectorOldOrder((1:numberTEs)+(iPLD-1)*numberTEs) = (iPLD-1)+1:numberPLDs:size(imASL,4);
                end
                imASL(:,:,:,1:end) = imASL(:,:,:,vectorOldOrder);
                xASL_io_SaveNifti(nii_files{1},nii_files{1},imASL);
                % Repeat Echo Times
                if numel(unique(resultJSON.EchoTime))>1 % Repeat EchoTime only if multiple TEs
                    resultJSON.EchoTime = repmat(resultJSON.EchoTime,numberPLDs,1);
                end
                % Save the JSON with the updated echo times
                spm_jsonwrite(fullfile(resultPath, [resultFile '.json']),resultJSON);
            else
                % Feedback about sequence
                warning('Time encoded sequence with one PLD only...');
            end
        end
    end

end

