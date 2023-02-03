function xASL_imp_DCM2NII_ReorderTimeEncoded(nii_files, bTimeEncoded, bTimeEncodedFME, timeEncodedMatrixSize, vectorPLD, resultJSON)
%xASL_imp_DCM2NII_ReorderTimeEncoded Reorder TEs and PLDs accordingly for time encoded sequences
%
% FORMAT: xASL_imp_DCM2NII_ReorderTimeEncoded(nii_files, bTimeEncoded, timeEncodedMatrixSize, vectorPLD, resultJSON)
% 
% INPUT:
%  nii_files     - List of nifti files (REQUIRED)
%  bTimeEncoded  - Time encoded sequence or not (BOOLEAN, REQUIRED)
%  bTimeEncoded  - Time encoded sequence from FME or not (BOOLEAN, REQUIRED)
%  timeEncodedMatrixSize - Time encoded matrix size (INTEGER or EMPTY, REQUIRED)
%  vectorPLD     - PLD vector (INTEGER or EMPTY, REQUIRED)
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
        % Check if the current sequence is a Hadamard from FME or not
		if bTimeEncoded && bTimeEncodedFME
            % Check image
			if xASL_exist(nii_files{1},'file')
                % Determine the number of time points within each NIfTI
                imASL = xASL_io_Nifti2Im(nii_files{1});
                numberTEs = length(unique(resultJSON.EchoTime));
				
				interleavedPLDs = false;
				if ~isempty(timeEncodedMatrixSize)
					numberPLDs = timeEncodedMatrixSize;
					numberRepetitions = int32(size(imASL,4)/numberTEs/numberPLDs);
					% Check if the PLDs are interleaved or just repeated
					if numel(vectorPLD) > numberPLDs
						for iRepetition = 2:floor(numel(vectorPLD)/numberPLDs)
							if ~isequal(vectorPLD(1:numberPLDs),vectorPLD((1:numberPLDs)+numberPLDs*(iRepetition-1)))
								interleavedPLDs = true;
							end
						end
					end
				else
					numberPLDs = int32(size(imASL,4)/numberTEs);
					numberRepetitions = 1;
				end
                
				if numberTEs > 1 && numberRepetitions > 0
					if numberRepetitions > 1 && interleavedPLDs
						error('Import of FME TimeEncoded for multiple TEs and Repetitions is not yet implemented for interleaved PLDs');
					end
					
					% Reorder TEs and PLDs - first cycle TE afterwards PLD
					vectorOldOrder = zeros(size(imASL,4),1);
					for iPLD = 1:(double(numberPLDs*numberRepetitions))
						vectorOldOrder((1:numberTEs)+(iPLD-1)*numberTEs) = (iPLD-1)+1:(numberPLDs*numberRepetitions):size(imASL,4);
					end
					imASL(:,:,:,1:end) = imASL(:,:,:,vectorOldOrder);
					xASL_io_SaveNifti(nii_files{1}, nii_files{1}, imASL);
					
					% Repeat Echo Times
					resultJSON.EchoTime = repmat(resultJSON.EchoTime, (numberPLDs*numberRepetitions), 1);
					
					% Save the JSON with the updated echo times
					xASL_io_WriteJson(fullfile(resultPath, [resultFile '.json']),resultJSON);
				elseif numberRepetitions > 1
					if ~interleavedPLDs 
						error('Import of FME TimeEncoded with single-TE and multi-Repetitions not yet implemented for non-interleaved PLDs');
					end
					% Reorder Repetitions and PLDs - first cycle PLDs afterwards Repetitions
					vectorOldOrder = zeros(size(imASL,4),1);
					for iRepetition = 1:(double(numberRepetitions))
						vectorOldOrder((1:numberPLDs)+(iRepetition-1)*numberPLDs) = (iRepetition-1)+1:numberRepetitions:size(imASL,4);
					end
					imASL(:,:,:,1:end) = imASL(:,:,:,vectorOldOrder);
					xASL_io_SaveNifti(nii_files{1},nii_files{1},imASL);
				end
                
            else
                % Feedback about sequence
                warning('Time encoded sequence with one PLD only...');
            end
        end
    end

end

