function xASL_adm_AlignLesionROItoStructural(Lesion_ROI_list)
%xASL_adm_AlignLesionROItoStructural Go through the list of Lesions and ROIs and double-check if they are aligned with their corresponding structural file
%
% FORMAT: xASL_adm_AlignLesionROItoStructural(Lesion_ROI_list)
%
% INPUT:
%   Lesion_ROI_list - list of Lesions and ROIs (REQUIRED)
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Go through all the Lesions and ROIs in the list. Check the corresponding FLAIR/T1 files. If the MAT of both files are similar, then do not do anything.
% if they differ, but the matrix size and MAT0 are equal, then reset the MAT for the Lesion or ROI and save again
%
% EXAMPLE: xASL_adm_AlignLesionROItoStructural(Lesion_ROI_list)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

if nargin<1
	error('Need 1 input argument');
end

% Cycle through all files
for iS=1:length(Lesion_ROI_list)
	% Obtain the structural file and check its existence
	[fPath, fName, fExt] = xASL_fileparts(Lesion_ROI_list{iS}); % Split the filename to path and name
	[startIndex,endIndex] = regexp(fName, '(T1|FLAIR)'); % Extract the name of the structural file

    correctIndexStart = ~isempty(startIndex) && isnumeric(startIndex) && all(isfinite(startIndex)) && numel(startIndex)==1;
    correctIndexEnd = ~isempty(endIndex) && isnumeric(endIndex) && all(isfinite(endIndex)) && numel(endIndex)==1;
    
    if correctIndexStart && correctIndexEnd

	    fStructName = fullfile(fPath,[fName(startIndex:endIndex) fExt]);
    
	    if xASL_exist(fStructName, 'file')
		    % Load NIfTI header of both files
		    lesionHeader = xASL_io_ReadNifti(Lesion_ROI_list{iS});
		    structHeader = xASL_io_ReadNifti(fStructName);
    
		    % Check if the MAT are OK
		    if ~all(abs(lesionHeader.mat - structHeader.mat) < 1e-3, 'all')
			    % The transformation matrices differ, we have to fix this
			    
			    if ~isequal(size(lesionHeader.dat), size(structHeader.dat))
				    % Option 2 - image sizes differ, we report a difference that we cannot fix
				    warning('%s\n%s\n%s', 'The transformation matrix MAT and the image size of Lesion/ROI and the corresponding structural files differ. Please check the alignment of the files:', Lesion_ROI_list{iS}, fStructName);
			    elseif ~all(abs(lesionHeader.mat0 - structHeader.mat0) < 1e-3, 'all')
				    % Option 3 - MAT0 also differ, we report a difference that we cannot fix
				    warning('%s\n%s\n%s','The transformation matrix MAT and MAT0 of Lesion/ROI and the corresponding structural files differ. Please check the alignment of the files:', Lesion_ROI_list{iS}, fStructName);
			    else
				    % Option 4 - we report a difference and set MAT of the Lesion/ROI to that of T1
				    fprintf('%s\n%s\n%s\n%s\n%s','The transformation matrix MAT of Lesion/ROI:', Lesion_ROI_list{iS}, ...
					                       'and the corresponding structural file:', fStructName, ...
									       ['differ, but MAT0 are equal. This means that the alignment is currently wrong, but was initially correct. The reason could be that Lesion/ROI was added after ExploreASL ' ...
									       'was executed. We have now set the Lesion/ROI orientation (==MAT) to that of the structural file. Please ensure that the files are correctly aligned.']);
    
				    % Save the Lesion/ROI again with a correct MAT and assign a warning
				    imLesion = xASL_io_Nifti2Im(Lesion_ROI_list{iS});
				    xASL_io_SaveNifti(fStructName, Lesion_ROI_list{iS}, imLesion);
    
				    error('Misalignment detected, tried to fix this, but a user-check is required before rerunning! See details above.');
			    end
            end
        end
    end
end

end

