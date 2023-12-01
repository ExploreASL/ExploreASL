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
% 2. Reorder multi-TE data
% 3. Decode the Hadamard data
% 4. Reorder multi-TE back to the initial order of PLD/TE
% 5. Normalization of the decoded data
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

%% 1. Specify the Encoding matrix

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
if length(xQ.LabelingDuration)>1 && length(xQ.LabelingDuration) ~= nEncodedVolumes
	error('Length of LabelingDuration vector does not match the number of volumes.');
end

if length(xQ.Initial_PLD)>1 && length(xQ.Initial_PLD) ~= nEncodedVolumes
	error('Length of Initial_PLD vector does not match the number of volumes.');
end

% Prepare the parameters also for interleaved reordering
nRepetitions = nEncodedVolumes / (xQ.TimeEncodedMatrixSize * xQ.NumberEchoTimes); % Calculating no. of acquisition repeats accounting for multipleTE

if length(xQ.EchoTime)>1 && length(xQ.EchoTime) ~= nEncodedVolumes
	error('Length of EchoTime vector does not match the number of volumes.');
end

% Currently, we only implement the option that TE times are together
% In case a different ordering is used, an error is reported
if isfield(xQ, 'NumberEchoTimes') && (xQ.NumberEchoTimes > 1)
	% The first two EchoTimes are equal
	if length(uniquetol(xQ.EchoTime(1:2), 0.001)) == 1
		error('Processing of multi-TE sequence is only implemented for the case of TEs being together in ASL4D.')
	end
end

%% 3. Obtain control4D images
% Here we select all TEs and the control images
% We thus calculate the size of each Hadamard block as the number of Hadamard phases and TEs
nVolumesPerRepetition = xQ.TimeEncodedMatrixSize * xQ.NumberEchoTimes;

% For example for 64 volumes and 2 repetitions with 8 PLDs and 4 TEs, it takes volume 1,2,3,4 and 33,34,35,36 to get all TEs of the control
indexControl = ((1:xQ.NumberEchoTimes)'*ones(1,nRepetitions) + ones(xQ.NumberEchoTimes,1)*(0:nRepetitions-1)*nVolumesPerRepetition);
indexControl = indexControl(:)';

decodedControl4D = imEncoded(:,:,:, indexControl);
xQ.EchoTime_Control4D = xQ.EchoTime(indexControl);
xQ.InitialPLD_Control4D = xQ.Initial_PLD(indexControl);
xQ.LabelingDuration_Control4D = xQ.LabelingDuration(indexControl);

%% 4. Reorder image matrix
% At this point the data is organized like this (in terms of ASL4D.nii volumes):
% PLD1/TE1,PLD1/TE2,PLD1/TE3,PLD1/TE4...PLD2/TE1,PLD2/TE2,PLD2/TE3... (TEs in first dimension, PLDs after)
%
% And for decoding we want
% TE1/PLD1,TE1/PLD2,TE1/PLD3,TE1/PLD4...TE2/PLD1,TE2/PLD2,TE2/PLD3,TE2/PLD4 (PLDs in the first dimension, TEs after)
% Then below, we apply this "transformation"/"decoding" vector to the image matrix and to the parameters TE/PLD/LD
if isfield(xQ, 'NumberEchoTimes') && (xQ.NumberEchoTimes > 1)
	% This is the original order
	vectorOldOrder = 1:nEncodedVolumes;
	
	% Shape to a matrix with TEs first and all the rest later
	vectorOldOrder = reshape(vectorOldOrder, xQ.NumberEchoTimes, nEncodedVolumes/xQ.NumberEchoTimes);
	
	% Flip the two dimensions and make a row vector again
	vectorOldOrder = reshape(vectorOldOrder', 1, nEncodedVolumes);

	% Reorder the data
	imEncoded = imEncoded(:,:,:,vectorOldOrder);
end

% Sometimes HAD8 is PLD1(rep1), PLD1(rep2), PLD2(rep1), PLD2(rep2). We
% want: PLD1(rep1), PLD2(rep1)... PLD1(rep2),PLD2(rep2)
% UniquePLDs = unique(xQ.Initial_PLD);
% 
% if (length(UniquePLDs) == nEncodedVolumes) && (UniquePLDs(2)-UniquePLDs(1))< 250
%     This is the original order
% 	vectorOldOrder = 1:nEncodedVolumes;
%     Shape
%     vectorOrder = [ vectorOldOrder(1:2:end) vectorOldOrder(2:2:end)];
%     Reorder the data
% 	imEncoded = imEncoded(:,:,:,vectorOrder);
% end

%% 5. Decode the Hadamard data
% We assume that the entire encoded volume is reorder into several repetitions of basic Hadamard blocks (repetitions or TE are counted the same)
decodedPWI4D = zeros(size(imEncoded,1), size(imEncoded,2), size(imEncoded,3), nEncodedVolumes / xQ.TimeEncodedMatrixSize * (xQ.TimeEncodedMatrixSize-1)); 
decodedBlock = zeros(size(imEncoded,1), size(imEncoded,2), size(imEncoded,3), xQ.TimeEncodedMatrixSize-1); % Temporary single decoded Hadamard block

% Go through all repetitions of the Hadamard blocks
for iRepetition = 1:(nEncodedVolumes / xQ.TimeEncodedMatrixSize)
	offsetEncoded = xQ.TimeEncodedMatrixSize * (iRepetition-1); % Index offset for the given repetition to search in the encoded data
	encodedBlock = imEncoded(:,:,:,(1:xQ.TimeEncodedMatrixSize) + offsetEncoded); % Obtain the single encoded block
	
	% Decoded the single block
	for iTI = 1:(xQ.TimeEncodedMatrixSize-1) %Number of decoded PLDs is number of Hadamard matrix size - 1
		indexPositive = find(xQ.TimeEncodedMatrix(iTI+1,:)==1);
		indexNegative = find(xQ.TimeEncodedMatrix(iTI+1,:)==-1);
		decodedBlock(:,:,:,iTI) = sum(encodedBlock(:,:,:,(indexPositive)),4) - sum(encodedBlock(:,:,:,(indexNegative)),4);
	end
	offsetDecoded = (iRepetition-1)*(xQ.TimeEncodedMatrixSize-1); % Index offset for the given repetition in the decoded data
	decodedPWI4D(:,:,:,offsetDecoded + (1:(xQ.TimeEncodedMatrixSize-1))) = decodedBlock; % Save the decoded block
end

%% 6. Reorder multi-TE back to the initial order of PLD/TE
% For model fitting, we want the PLDs-first-TEs-second order (just like at the beginning) so we need to reorder it again
if isfield(xQ,'NumberEchoTimes') && (xQ.NumberEchoTimes > 1)
	% This is the original order
	vectorOldOrder = 1:size(decodedPWI4D, 4);
	
	% Shape to a matrix with TEs second and all the rest first
	vectorOldOrder = reshape(vectorOldOrder, size(decodedPWI4D, 4)/xQ.NumberEchoTimes, xQ.NumberEchoTimes);
	
	% Flip the two dimensions and make a row vector again
	vectorOldOrder = reshape(vectorOldOrder', 1, size(decodedPWI4D, 4));

	% Reorder the data
	decodedPWI4D = decodedPWI4D(:,:,:,vectorOldOrder);
end

% Calculate the correct parameters
indexDecoding = ((xQ.NumberEchoTimes+1):nVolumesPerRepetition)'*ones(1,nRepetitions) + ones(nVolumesPerRepetition-xQ.NumberEchoTimes,1)*(0:nRepetitions-1)*nVolumesPerRepetition;
indexDecoding = indexDecoding(:)';
xQ.EchoTime_PWI4D = xQ.EchoTime(indexDecoding);
xQ.InitialPLD_PWI4D = xQ.Initial_PLD(indexDecoding);
xQ.LabelingDuration_PWI4D = xQ.LabelingDuration(indexDecoding);

%% 7. Normalization of the decoded data
% NormalizationFactor = 1/(m_INumSets/2);
% where m_INumSets is the number of images (e.g. 8 for Hadamard 8x8)

decodedPWI4D = decodedPWI4D / (xQ.TimeEncodedMatrixSize/2);
    
end