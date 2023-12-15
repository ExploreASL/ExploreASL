function [decodedPWI4D, decodedControl4D, xQ] = xASL_quant_HadamardDecoding(imPath, xQ)
%xASL_quant_HadamardDecoding Hadamard-4 & Hadamard-8 Decoding
%
% FORMAT:       [imDecoded] = xASL_quant_HadamardDecoding(imPath, xQ)
%
% INPUT:        imPath - Path to the encoded ASL4D image we want to decode (STRING OR 4D MATRIX, REQUIRED)
%                        Alternatively, this can contain the image matrix
%                        (xASL_io_Nifti2Im below allows both path and image
%                        matrix inputs)
%
%               xQ     - xQ field with Hadamard input parameters containing the following subfields
%                          - TimeEncodedMatrixType (REQUIRED)
%                             - Hadamard
%                             - Walsh
%                          - TimeEncodedMatrixSize (REQUIRED)
%                             - '4' for Hadamard-4
%                             - '8' for Hadamard-8
%                          - TimeEncodedMatrix (OPTIONAL)
%                             - Matrix given by the user
%                           - EchoTime: TE vector for encoded images
%                           - Initial_PLD: PLD vector for encoded images
%                           - LabelingDuration: : LD vector for encoded images
%
% OUTPUT:       decodedPWI4D - Decoded ASL volumes
%               decodedControl4D - Decoded control volumes
%               xQ     - xQ field with several output parameters
%                           - EchoTime_PWI4D: TE vector after Hadamard-decoded subtraction
%                           - InitialPLD_PWI4D: PLD vector after Hadamard-decoded subtraction
%                           - LabelingDuration_PWI4D: LD vector after Hadamard-decoded subtraction
%                           - EchoTime_Control4D: TE vector for all control images
%                           - InitialPLD_Control4D: PLD vector for all control images
%                           - LabelingDuration_Control4D: : LD vector for all control images
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Hadamard-4 & Hadamard-8 Decoding.
%
% 0. Admin: Check inputs, load data
% 1. Specify the decoding matrix
% 2. Load time-series nifti
% 3. Calculate the number of repetitions
% 4. Obtain control4D images
% 5. Reorder multi-TE data
% 6. Decode the image volumes according to the Hadamard scheme
% 7. Normalization of the decoded data
% 8. Reorder multi-TE back to the initial order of PLD/TE
% 9. Calculate the correct parameters of the decoded data
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      n/a
% __________________________________
% Copyright (c) 2015-2023 ExploreASL

%% 0. Admin
% Check if all inputs are present

if nargin<1 || isempty(imPath)
    error('imPath input is empty');
end

if nargin<2 || isempty(xQ)
    error('xQ input is empty');
end

%% 1. Specify the decoding matrix

% Specify the TimeEncodedMatrix
if ~isempty(xQ.TimeEncodedMatrix)
	if size(xQ.TimeEncodedMatrix,1) ~= size(xQ.TimeEncodedMatrix,2)
		error('TimeEncodedMatrix must be square NxN ');
	end
	
    % TimeEncodedMatrix exists, we verify the size
	if isempty(xQ.TimeEncodedMatrixSize)
		xQ.TimeEncodedMatrixSize = size(xQ.TimeEncodedMatrix,2);
	else
		if size(xQ.TimeEncodedMatrix,2) ~= xQ.TimeEncodedMatrixSize
			error('Mismatch between TimeEncodedMatrix and TimeEncodedMatrixSize');
		end
	end
end

if ~isempty(x.Q.TimeEncodedMatrixType)
	% See an example of decoding/encoding matrices in 
	% Samson-Himmelstjerna, MRM 2015 https://doi.org/10.1002/mrm.26078
	% Note that the encoding/decoding matrices are symmetric
	% Encoding matrices have an extra first row with all 1 (1 control, -1 label)
	% Decoding can be done using the same matrix (1 addition, -1 subtraction) - below, the first row 
	% (that would generate a mean control) is skipped and then we use the decoding along columns, goind across
	% rows generates different decoded volumes.
	% #### For Walsh Decoding Matrix ####
	if strcmpi(xQ.TimeEncodedMatrixType,'Walsh')
		
		if xQ.TimeEncodedMatrixSize == 4
			tempTimeEncodedMatrix =...
				[1  1  1  1; 
				 1 -1  1 -1;
				 1 -1 -1  1;
				 1  1 -1 -1];
		elseif xQ.TimeEncodedMatrixSize == 8
            tempTimeEncodedMatrix = [1  1  1  1  1  1  1  1;
				                     1 -1  1 -1  1 -1  1 -1;
                                     1 -1  1 -1 -1  1 -1  1;
                                     1 -1 -1  1 -1  1  1 -1;
                                     1 -1 -1  1  1 -1 -1  1;
                                     1  1 -1 -1  1  1 -1 -1;
                                     1  1 -1 -1 -1 -1  1  1;
                                     1  1  1  1 -1 -1 -1 -1];
		else
			tempTimeEncodedMatrix = [];
		end
		
		% #### For Hadamard Decoding Matrix ####
	elseif strcmpi(xQ.TimeEncodedMatrixType,'Hadamard')
		% An alternative Philips version xQ.TimeEncodedMatrix 
		% HAD4 [1, 1, 1, 1;
		%       1, 1,-1,-1;
		%       1,-1, 1,-1;
		%       1,-1,-1,1];
		% HAD8 [1, 1, 1, 1, 1, 1, 1, 1;
		%       1,-1, 1,-1, 1,-1, 1,-1;
		%       1, 1,-1,-1, 1, 1,-1,-1;
		%       1,-1,-1, 1, 1,-1,-1, 1;
		%       1, 1, 1, 1,-1,-1,-1,-1;
		%       1,-1, 1,-1,-1, 1,-1, 1;
		%       1, 1,-1,-1,-1,-1, 1, 1;
		%       1,-1,-1, 1,-1, 1, 1,-1];
		if xQ.TimeEncodedMatrixSize == 4
            tempTimeEncodedMatrix = [1  1  1  1;
				                     1 -1  1 -1;
                                     1  1 -1 -1;
                                     1 -1 -1  1];
		elseif xQ.TimeEncodedMatrixSize == 8
            tempTimeEncodedMatrix = [1  1  1  1  1  1  1  1;
				                     1 -1 -1  1 -1  1  1 -1;
                                     1  1 -1 -1 -1 -1  1  1;
                                     1 -1  1 -1 -1  1 -1  1;
                                     1  1  1  1 -1 -1 -1 -1;
                                     1 -1 -1  1  1 -1 -1  1;
                                     1  1 -1 -1  1  1 -1 -1;
                                     1 -1  1 -1  1 -1  1 -1];
		else
			tempTimeEncodedMatrix = [];
		end    
	end
	if isempty(tempTimeEncodedMatrix)
		warning(['Cannot create a ' xQ.TimeEncodedMatrixType ' matrix of size ' xASL_num2str(xQ.TimeEncodedMatrixSize)]);
	end

	if isempty(xQ.TimeEncodedMatrix)
		% If the matrix is not yet initialize, then save it
		xQ.TimeEncodedMatrix = tempTimeEncodedMatrix;
	else
		% If the created and initialized matrices are equal
		if ~isequal(tempTimeEncodedMatrix, xQ.TimeEncodedMatrix)
			% If the created and initialized matrices differ, then issue a warning and use the user matrix
			warning('Created matrix based on provided Type and size differs from the provided Matrix. Using the provided matrix');
		end
	end
end
   
switch xQ.Vendor
	case 'Siemens'
		% This is where this code was initially based on, so for
		% Siemens we don't do anything specific
	case 'Philips'
		xQ.TimeEncodedMatrix = xQ.TimeEncodedMatrix * -1;
		% Philips uses an inverse matrix compared to Siemens
	case 'GE'
		warning('Time encoded ASL data detected for GE, but GE usually does the decoding in K-space on the scanner');
	otherwise
		warning(['Time encoded ASL data not yet tested for vendor ' xQ.vendor]);
end


%% 2. Load time-series nifti
imEncoded = xASL_io_Nifti2Im(imPath); 

% Prepare the general sequence parameters
nEncodedVolumes = size(imEncoded, 4);     % Size of the raw data

% Check that input quantification parameters are vectors with an equal length to the number of volumes
if length(xQ.LabelingDuration) ~= nEncodedVolumes
	if length(xQ.LabelingDuration) == 1
		xQ.LabelingDuration = ones(1, nEncodedVolumes) * xQ.LabelingDuration;
	else
		error('Length of LabelingDuration vector does not match the number of volumes.');
	end
end

if length(xQ.Initial_PLD) ~= nEncodedVolumes
	if length(xQ.Initial_PLD) == 1
		xQ.Initial_PLD = ones(1, nEncodedVolumes) * xQ.Initial_PLD;
	else
		error('Length of Initial_PLD vector does not match the number of volumes.');
	end
end

if length(xQ.EchoTime) ~= nEncodedVolumes
	if length(xQ.EchoTime) == 1
		xQ.EchoTime = ones(1, nEncodedVolumes) * xQ.EchoTime;
	elseif mod(nEncodedVolumes, length(xQ.EchoTime)) == 0
		xQ.EchoTime = repmat(xQ.EchoTime(:), [nEncodedVolumes/length(xQ.EchoTime) 1]);
	else
		error('Length of EchoTime vector does not match the number of volumes.');
	end
end

% Currently, for multi-TE option, we only implement the option that TE times are in the first dimension
% So in case that there's multi-TE but not as the inner most dimension, we report an error
if isfield(xQ, 'NumberEchoTimes') && (xQ.NumberEchoTimes > 1)
	% The first two EchoTimes are equal
	if length(uniquetol(xQ.EchoTime(1:2), 0.001)) == 1
		warning('Processing of multi-TE sequence is only implemented for the case of TEs being the inner-most dimension in ASL4D. The current option might not work correctly.')
	end
end

%% 3. Calculate the number of repetitions
% Here we define the terminology about repetitions and blocks
% 1. We have a set of several individual volumes that we call nEncodedVolumes
% 2. Several volumes - i.e. several blocks - are acquired together corresponding to a fully encoded Hadamard matrix (e.g. 4 or 8 blocks together)
%    we call the number of the nFullHadamardVolumes
nFullHadamardVolumes = nEncodedVolumes / xQ.TimeEncodedMatrixSize;
% 3. ASL measurements are typically repeated several times. Though note that with Hadamard encoded we acquire several volumes for a full Hadamard volume, and for segmented 3D readout, the full 3D
%    volume is segmented into several readouts. So for 3D Grase Hadamard, there's typically only time to make 1-3 full repetitions. Besides full repetitions that refer to repeating the entire acquisition scheme,
%    we acquire the Hadamard volume multiple times in a row for other reasons. Typically, a multi-echo acquisition is done. This results in several acquired Hadamard volumes without really repeating the whole 
%    acquisition process. We call this here "inner repetitions". Normally, the PLD changes for each Hadamard block, so if there are several equal PLDs together, we identify this as an inner repetition.
%    This could be inner-ly repeated PLDs for some implementation reason or because of multiTE. We check both and issue a warning if more PLDs than TEs are present, as we don't really implement this option
%    properly now. 
[~,~,indexFirstPLD] = unique(xQ.Initial_PLD); % Find locations of all PLDs
indexSecondPLD = find(indexFirstPLD==2); % Find the location of the second PLD
if isempty(indexSecondPLT)
	nInnerRepetitions = length(xQ.Initial_PLD); % There's no second PLD and nInnerRepetitions is the length of the PLD vector
else
	nInnerRepetitions = indexSecondPLD - 1; % The number of first-PLDs, i.e. the count of all PLDs before the second one, is the number of inner repetitions
end

if nInnerRepetitions ~= xQ.NumberEchoTimes % Check for TE number differ from the number of inner repeats as this is a special option needing a special reordering
	warning('Number of inner repetitions and echo-times does not match for multi-TE Hadamard. This option is not yet properly implemented.');
end

% 4. We define how many volumes are within an inner repetition accounting for the number of blocks within a Hadamard volume
nVolumesPerInnerRepetition = xQ.TimeEncodedMatrixSize * nInnerRepetitions;

% 5. We calculate the number of full repetitions - that is reapeating the entire acquisition scheme
nFullRepetitions = nEncodedVolumes / nVolumesPerInnerRepetition; 


%% 4. Obtain control4D images
% Here we select all TEs for the control images
% For example for 64 volumes and 2 repetitions with 8 PLDs and 4 TEs, it takes volume 1,2,3,4 and 33,34,35,36 to get all TEs of the control
indexControl = ((1:nInnerRepetitions)'*ones(1,nFullRepetitions) + ones(nInnerRepetitions,1)*(0:nFullRepetitions-1)*nVolumesPerInnerRepetition);
indexControl = indexControl(:)';

decodedControl4D = imEncoded(:,:,:, indexControl);
xQ.EchoTime_Control4D = xQ.EchoTime(indexControl);
xQ.InitialPLD_Control4D = xQ.Initial_PLD(indexControl);
xQ.LabelingDuration_Control4D = xQ.LabelingDuration(indexControl);

%% 5. Reorder image matrix
% Reorders that data for innerRepetitions. The example below shows it for the case of innerRepetitions being TEs only:
% At this point the data is organized like this (in terms of ASL4D.nii volumes):
% PLD1/TE1,PLD1/TE2,PLD1/TE3,PLD1/TE4...PLD2/TE1,PLD2/TE2,PLD2/TE3... (TEs in first dimension, PLDs after)
%
% And for decoding we want to swap the dimensions for TE and PLD so that:
% TE1/PLD1,TE1/PLD2,TE1/PLD3,TE1/PLD4...TE2/PLD1,TE2/PLD2,TE2/PLD3,TE2/PLD4 (PLDs in the first dimension, TEs after)
% Then below, we apply this "transformation"/"decoding" vector to the image matrix and to the parameters TE/PLD/LD
if nInnerRepetitions > 1)
	% This is the original order
	vectorReorderEncoded = 1:nEncodedVolumes;
	
	% Shape to a matrix with TEs first and all the rest later
	vectorReorderEncoded = reshape(vectorReorderEncoded, nInnerRepetitions, nEncodedVolumes/nInnerRepetitions);
	
	% Switch the two dimensions and create a row vector again
	vectorReorderEncoded = reshape(vectorReorderEncoded', 1, nEncodedVolumes);

	% Reorder the data
	imEncodedReordered = imEncoded(:,:,:,vectorReorderEncoded);
end

% Sometimes HAD8 is PLD1(rep1), PLD1(rep2), PLD2(rep1), PLD2(rep2). We
% want: PLD1(rep1), PLD2(rep1)... PLD1(rep2),PLD2(rep2)

%% 6. Decode the image volumes according to the Hadamard scheme
% We assume that the entire encoded volume is reorder into several repetitions of basic Hadamard blocks
% This is not a true repetition in the sense of repeating the entire acquisition as it counts repeats also for different TEs.
% So we call it nFullHadamardVolumes

decodedPWI4D = zeros(size(imEncodedReordered,1), size(imEncodedReordered,2), size(imEncodedReordered,3), nFullHadamardVolumes * (xQ.TimeEncodedMatrixSize-1)); 

% Go through all repetitions of the full Hadamard blocks
for iRepetition = 1:nFullHadamardVolumes
	offsetEncoded = xQ.TimeEncodedMatrixSize * (iRepetition-1); % Index offset for the given repetition to search in the encoded data
	encodedBlocks = imEncoded(:,:,:,(1:xQ.TimeEncodedMatrixSize) + offsetEncoded); % Obtain the single encoded block
	decodedBlocks = zeros(size(imEncodedReordered,1), size(imEncodedReordered,2), size(imEncodedReordered,3), xQ.TimeEncodedMatrixSize-1); % Temporary single decoded Hadamard block
	
	% Decoded the single block
	for iBlock = 1:(xQ.TimeEncodedMatrixSize-1) %Number of decoded PLDs is number of Hadamard matrix size - 1
		indexPositive = find(xQ.TimeEncodedMatrix(iBlock+1,:)==1);
		indexNegative = find(xQ.TimeEncodedMatrix(iBlock+1,:)==-1);
		decodedBlocks(:,:,:,iBlock) = sum(encodedBlocks(:,:,:,(indexPositive)),4) - sum(encodedBlocks(:,:,:,(indexNegative)),4);
	end
	offsetDecoded = (iRepetition-1)*(xQ.TimeEncodedMatrixSize-1); % Index offset for the given repetition in the decoded data
	decodedPWI4D(:,:,:,offsetDecoded + (1:(xQ.TimeEncodedMatrixSize-1))) = decodedBlocks; % Save the decoded block
end

%% 7. Normalization of the decoded data
% NormalizationFactor = 1/(m_INumSets/2);
% where m_INumSets is the number of images (e.g. 8 for Hadamard 8x8)

decodedPWI4D = decodedPWI4D / (xQ.TimeEncodedMatrixSize/2);

%% 8. Reorder multi-TE back to the initial order of PLD/TE
% For model fitting, we want the PLDs-first-TEs-second order (just like at the beginning) so we need to reorder it again
if nInnerRepetitions > 1
	% This is the original order
	vectorReorderDecoded = 1:size(decodedPWI4D, 4);
	
	% Shape to a matrix with TEs second and all the rest first
	vectorReorderDecoded = reshape(vectorReorderDecoded, size(decodedPWI4D, 4)/nInnerRepetitions, nInnerRepetitions);
	
	% Flip the two dimensions and make a row vector again
	vectorReorderDecoded = reshape(vectorReorderDecoded', 1, size(decodedPWI4D, 4));

	% Reorder the data
	decodedPWI4D = decodedPWI4D(:,:,:,vectorReorderDecoded);
end

%% 9. Calculate the correct parameters of the decoded data
% Calculate the correct parameters
indexDecoding = ((nInnerRepetitions+1):nVolumesPerInnerRepetition)'*ones(1,nFullRepetitions) + ones(nVolumesPerInnerRepetition-nInnerRepetitions,1)*(0:nFullRepetitions-1)*nVolumesPerInnerRepetition;
indexDecoding = indexDecoding(:)';
xQ.EchoTime_PWI4D = xQ.EchoTime(indexDecoding);
xQ.InitialPLD_PWI4D = xQ.Initial_PLD(indexDecoding);
xQ.LabelingDuration_PWI4D = xQ.LabelingDuration(indexDecoding);
    
end