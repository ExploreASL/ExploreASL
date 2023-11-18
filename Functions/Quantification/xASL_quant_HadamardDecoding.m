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
%
%                           - EchoTime: TE vector after Hadamard-decoded subtraction
%                           - Initial_PLD: PLD vector after Hadamard-decoded subtraction
%                           - LabelingDuration: LD vector after Hadamard-decoded subtraction
%                           - EchoTime: TE vector for all control images
%                           - Initial_PLD: PLD vector for all control images
%                           - LabelingDuration: : LD vector for all control images
%
% OUTPUT:       decodedPWI4D - Decoded ASL volumes
%               decodedControl4D - Decoded control volumes
%               xQ     - xQ field with several output parameters
%                           - HadamardDecodingVector (REQUIRED)
%                             vector with element-wise index where a volume should go
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

% Check that quantification parameters are vectors with an equal length to the number of volumes
% otherwise we cannot decode them


% if length(xQ.EchoTime)

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

% Calculate specific data sizes
nDecodedTI = xQ.TimeEncodedMatrixSize-1;                   % Number of TIs is always MatrixSize -1

%% #997 Here we need to use something different than xQ.NumberEchoTimes to determine the structure
%  Time decoding should be oblivious about the readout parameters
nDecodedVolumes = xQ.NumberEchoTimes * nDecodedTI;

nEncodedVolumes = size(imEncoded, 4);

% #997 Here we need to use something different than xQ.NumberEchoTimes to determine the structure
nRepetitions = nEncodedVolumes / (xQ.TimeEncodedMatrixSize * xQ.NumberEchoTimes); % Calculating no. of acquisition repeats

%% 3. Obtain control4D images
% Here we select the 1st TE and the control images
% We thus calculate the size of each Hadamard block as the number of Hadamard phases and TEs
nVolumesPerRepetition = xQ.TimeEncodedMatrixSize * xQ.NumberEchoTimes;
% #997 Here we need to use something different than xQ.NumberEchoTimes to determine the structure

% For example for 64 volumes and 2 repetitions with 8 PLDs and 4 TEs, it takes volume 1 and 33
decodedControl4D = imEncoded(:,:,:,1:nVolumesPerRepetition:end);

%% #997 Here, we need to create a 1D decoding vector, that for each element (original volume) has the new position of the volume
%  Volumes that disappear will have a 0
% e.g., 3 3 3 3 2 2 2 2 1 1 1 1 0 0 0 0
% Then below, we apply this "transformation"/"decoding" vector to the image matrix and to the parameters TE/PLD/LD


%% 4. Reorder image matrix
% At this point the data is organized like this (in terms of ASL4D.nii volumes):
% PLD1/TE1,PLD1/TE2,PLD1/TE3,PLD1/TE4...PLD2/TE1,PLD2/TE2,PLD2/TE3... (TEs in first dimension, PLDs after)
%
% And for decoding we want
% TE1/PLD1,TE1/PLD2,TE1/PLD3,TE1/PLD4...TE2/PLD1,TE2/PLD2,TE2/PLD3,TE2/PLD4 (PLDs in the first dimension, TEs after)
if isfield(xQ,'NumberEchoTimes') && (xQ.NumberEchoTimes > 1)
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
decodedPWI4D = zeros(size(imEncoded,1), size(imEncoded,2), size(imEncoded,3), nDecodedVolumes * nRepetitions);    
idx=0;
for iRepetition = 1:nRepetitions
    for iTE = 1:xQ.NumberEchoTimes
        for iTI = 1:nDecodedTI
            indexPositive = find(xQ.TimeEncodedMatrix(iTI+1,:)==1);
            indexNegative = find(xQ.TimeEncodedMatrix(iTI+1,:)==-1);
			decodedPWI4D(:,:,:,((iTE-1)*nDecodedTI+iTI)+(iRepetition-1)*nDecodedVolumes) = sum(imEncoded(:,:,:,(indexPositive+idx)),4) - sum(imEncoded(:,:,:,(indexNegative+idx)),4);
        end
        idx = idx+xQ.TimeEncodedMatrixSize;
    end
end


%% 6. Reorder multi-TE back to the initial order of PLD/TE
%
% For model fitting, we want the PLDs-first-TEs-second order (just like the
% beginning) so we need to reorder it again
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


%% 7. Normalization of the decoded data
% NormalizationFactor = 1/(m_INumSets/2);
% where m_INumSets is the number of images (e.g. 8 for Hadamard 8x8)

decodedPWI4D = decodedPWI4D / (xQ.TimeEncodedMatrixSize/2);
    
end