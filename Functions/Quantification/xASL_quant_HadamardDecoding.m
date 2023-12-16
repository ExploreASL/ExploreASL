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
% 2. Verify that quantification parameters are equally long as the number of volumes
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
else
    imEncoded = xASL_io_Nifti2Im(imPath);
    nEncodedVolumes = size(imEncoded, 4); % Number of volumes of the raw data
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
	% (which would generate control volumes) is skipped and then we process the decoding matrix column-wise,
	% because processing row-wise would generate different/wrongly decoded volumes.
    %
	% #### For data encoded according to a Walsh matrix ####
	if strcmp(xQ.TimeEncodedMatrixType, 'Walsh')
		
		if xQ.TimeEncodedMatrixSize == 4
			tempTimeEncodedMatrix =...
				[1  1  1  1; % control-label
 				 1 -1  1 -1; 
				 1 -1 -1  1;
				 1  1 -1 -1];

		elseif xQ.TimeEncodedMatrixSize == 8
            tempTimeEncodedMatrix = [1  1  1  1  1  1  1  1;% control-label
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
		
		% #### For data encoded according to a Hadamard matrix ####
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
            tempTimeEncodedMatrix = [1  1  1  1;% control-label
				                     1 -1  1 -1;
                                     1  1 -1 -1;
                                     1 -1 -1  1];
		elseif xQ.TimeEncodedMatrixSize == 8
            tempTimeEncodedMatrix = [1  1  1  1  1  1  1  1;% control-label
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
	case 'GE'
		warning('Time encoded ASL data detected for GE, but GE usually does the decoding in K-space on the scanner. So consider disabling Hadamard decoding in ExploreASL');
	case 'Philips'
		xQ.TimeEncodedMatrix = xQ.TimeEncodedMatrix * -1;
		% Philips uses an inverse matrix compared to the original mathematical definitions
	case 'Siemens'
		% Siemens uses the original mathematical definition of Hadamard & Walsh matrices
	otherwise
		error(['Time encoded ASL data not implemented yet for vendor ' xQ.vendor]);
end

%% 2. Verify that quantification parameters are equally long as the number of volumes

if length(xQ.LabelingDuration) ~= nEncodedVolumes
	if length(xQ.LabelingDuration) == 1
        % If a single LD is specified, we assume that all volumes have the same LD
		xQ.LabelingDuration = ones(1, nEncodedVolumes) * xQ.LabelingDuration;
	else
		error('Length of LabelingDuration vector does not match the number of volumes');
	end
end

if length(xQ.Initial_PLD) ~= nEncodedVolumes
	if length(xQ.Initial_PLD) == 1
        % If a single PLD is specified, we assume that all volumes have the same PLD
		xQ.Initial_PLD = ones(1, nEncodedVolumes) * xQ.Initial_PLD;
	else
		error('Length of Initial_PLD vector does not match the number of volumes');
	end
end

if length(xQ.EchoTime) ~= nEncodedVolumes
	if length(xQ.EchoTime) == 1
        % If a single TE is specified, we assume that all volumes have the same TE
		xQ.EchoTime = ones(1, nEncodedVolumes) * xQ.EchoTime;
	elseif mod(nEncodedVolumes, length(xQ.EchoTime)) == 0
		xQ.EchoTime = repmat(xQ.EchoTime(:), [nEncodedVolumes/length(xQ.EchoTime) 1]);
	else
		error('Length of EchoTime vector does not match the number of volumes');
	end
end

% Currently, for multi-TE, we only implement the option that TE times are in the first dimension
% So in case that there's multi-TE but not as the inner most dimension, we report an error
if isfield(xQ, 'NumberEchoTimes') && (xQ.NumberEchoTimes > 1)
	% The first two EchoTimes are equal
	if length(uniquetol(xQ.EchoTime(1:2), 0.001)) == 1
		warning('Multi-TE processing is only implemented TEs being the inner-most dimension only in PWI4D repetitions, hence the current implementation may not work correctly for the current data')
	end
end


%% 3. Calculate the number of repetitions

% #997 -> Jan I tried to rewrite it, so it becomes fully understandable.
% Why don't we actually use the names nBlocks and nVolumesPerBlock
% isn't that much less confusion-prone? Because "nFullHadamardVolumes" to me could mean different things

% Here we define the terminology about repetitions and blocks
% nEncodedVolumes = total original volumes
%
%   A "block" is in mathematics a submatrix, a section of a matrix
%   With Hadamard-encoded ASL, the fully encoded Hadamard matrix is separated in e.g., 4 or 8 blocks (corresponding to Hadamard-4 or Hadamard-8, with
%   1 control volume and 3 or 7 PLDs, respectively). So:
%
% block = control or single PLD acquisition encoded in Hadamard
% TimeEncodedMatrixSize = nBlocks

nBlocks = xQ.TimeEncodedMatrixSize;

%   INNER REPETITIONS
%   There may be readout-related reasons for having multiple volumes.
%   E.g., with a multi-TE readout, each control or label volume (in this case each Hadamard-encoded volume) is acquired multiple 
%   times for the number of echoes. Because this is repeated at the readout level, we define these here as "inner repetitions".
%   E.g., for a 8 echoes Hadamard acquisition, a single control or PLD-encoded acquisition is repeated 8 times, so:
%
% innerRepetition: single-volume repetitions within a Hadamard block

[~, ~, indexUniquePLD] = unique(xQ.Initial_PLD); % Find volume indices of unique PLDs
indexFirstPLD = find(indexUniquePLD==1); % Find the volumes of the first PLD
indexSecondPLD = find(indexUniquePLD==2); % Find the volumes of the second PLD
if isempty(indexSecondPLD)
	nInnerRepetitions = length(xQ.Initial_PLD); % There's no second PLD and nInnerRepetitions is the length of the PLD vector
else
	nInnerRepetitions = indexSecondPLD - 1; % The number of first-PLDs, i.e. the count of all PLDs before the second one, is the number of inner repetitions
end

if nInnerRepetitions ~= xQ.NumberEchoTimes % Check for TE number differ from the number of inner repeats as this is a special option needing a special reordering
	warning('Number of inner repetitions and echo-times does not match, such a sequence is not yet implemented');
end

nVolumesPerBlock = nEncodedVolumes / nBlocks;

% We define how many volumes are within an inner repetition accounting for the number of blocks within a Hadamard volume
nVolumesPerInnerRepetition = nBlocks * nInnerRepetitions;

% We calculate the number of full repetitions - that is repeating the entire acquisition scheme
nAcquisitionRepetitions = nEncodedVolumes / nVolumesPerInnerRepetition; 

%   OUTER REPETITIONS
%   The full Hadamard-encoded combination of ASL measurements can be repeated for averaging to increase SNR. 
%   Though note that some readouts like 3D GRASE are typically segmented to reduce the acquisition PSF in the Z-direction, 
%   increasing the acquisition of a single control or label volume with a factor of 2, 4, 8 (depending on the number of segments).
%   Hence, for 3D Grase Hadamard, there's typically only clinically scanning time to scan 1-3 repetitions, which we here define as "outer repetitions",
%   so:
%
% outerRepetition: multi-volume repetitions containing full Hadamard-encoded ASL blocks

nOuterRepetitions = nVolumesPerBlock / nInnerRepetitions;


%% 4. Decode control4D
% Here we select all control volumes
% If there are multiple TEs or other readout-related repetitions, it takes all of them.
% For example for 64 volumes and 2 repetitions with 8 PLDs and 4 TEs, it takes volume 1,2,3,4 and 33,34,35,36 to get all control volumes
indexControl = ((1:nInnerRepetitions)'*ones(1,nAcquisitionRepetitions) + ones(nInnerRepetitions,1)*(0:nAcquisitionRepetitions-1)*nVolumesPerInnerRepetition);
indexControl = indexControl(:)';

decodedControl4D = imEncoded(:,:,:, indexControl);
xQ.EchoTime_Control4D = xQ.EchoTime(indexControl);
xQ.InitialPLD_Control4D = xQ.Initial_PLD(indexControl);
xQ.LabelingDuration_Control4D = xQ.LabelingDuration(indexControl);


%% 5. Switch Hadamard block dimension and innerRepetition dimension
% Prepares the ASL volumes for Hadamard decoding, by 
% reordering volumes such that Hadamard blocks are the most inner dimension/repetition, and any innerRepetitions (such as multi-TE) 
% are one dimension layer outside. 
%
% The following example illustrates this for a multi-TE innerRepetition:
% At this point the image volumes are organized as:
% PLD1/TE1,PLD1/TE2,PLD1/TE3,PLD1/TE4...PLD2/TE1,PLD2/TE2,PLD2/TE3... (TEs in first dimension, Hadamard blocks after)
%
% And for decoding we want to swap the dimensions such that the Hadamard blocks (==PLDs) are the innermost dimension:
% TE1/PLD1,TE1/PLD2,TE1/PLD3,TE1/PLD4...TE2/PLD1,TE2/PLD2,TE2/PLD3,TE2/PLD4 (blocks in the first dimension, TEs after)

% This is the original order
vectorReorderEncoded = 1:nEncodedVolumes;

% Shape to a matrix with TEs first and all the rest later
vectorReorderEncoded = reshape(vectorReorderEncoded, nInnerRepetitions, nEncodedVolumes/nInnerRepetitions);

% Switch the two dimensions and create a row vector again
vectorReorderEncoded = reshape(vectorReorderEncoded', 1, nEncodedVolumes);

% Reorder the data
imEncodedReordered = imEncoded(:,:,:,vectorReorderEncoded);


% Sometimes HAD8 is PLD1(rep1), PLD1(rep2), PLD2(rep1), PLD2(rep2). We
% want: PLD1(rep1), PLD2(rep1)... PLD1(rep2),PLD2(rep2)


%% 6. Decode the image volumes
% We assume that the all encoded volumes are ordered as repetitions of Hadamard blocks

nPLDs = nBlocks-1;
% We only decode the PLD blocks here, we ignore the control volumes
% Hence the number of PLDs is the total number of Hadamard blocks (TimeEncodedMatrixSize - 1)

decodedPWI4D = zeros(size(imEncodedReordered,1), size(imEncodedReordered,2), size(imEncodedReordered,3), nOuterRepetitions * nPLDs);

% Go through all repetitions of the full Hadamard blocks
for iRepetition = 1:nOuterRepetitions
	offsetEncoded = nBlocks * (iRepetition-1); % Index offset for the given repetition to search in the encoded data
	encodedBlocks = imEncoded(:,:,:,(1:nBlocks) + offsetEncoded); % Obtain the single encoded block
	decodedBlocks = zeros(size(imEncodedReordered,1), size(imEncodedReordered,2), size(imEncodedReordered,3), nBlocks-1); % Temporary single decoded Hadamard block

	for iBlock = 1:nPLDs

		indexPositive = find(xQ.TimeEncodedMatrix(iBlock+1,:)==1);
		indexNegative = find(xQ.TimeEncodedMatrix(iBlock+1,:)==-1);
		decodedBlocks(:,:,:,iBlock) = sum(encodedBlocks(:,:,:,(indexPositive)),4) - sum(encodedBlocks(:,:,:,(indexNegative)),4);
	end
	offsetDecoded = (iRepetition-1)*nPLDs; % Volume index offset for the current outer repetition
	decodedPWI4D(:,:,:,offsetDecoded + (1:nPLDs)) = decodedBlocks; % Save the decoded block
end

% Normalization of the decoded data
% NormalizationFactor = 1/(m_INumSets/2);
% where m_INumSets is the number of images (e.g. 8 for Hadamard 8x8)

decodedPWI4D = decodedPWI4D / (xQ.TimeEncodedMatrixSize/2);

%% 7. Switch Hadamard block dimension and innerRepetition dimension back
% For model fitting, we want the PLDs-first-TEs-second order (just like at the beginning) so we need to reorder it again

% This is the original order
vectorReorderDecoded = 1:size(decodedPWI4D, 4);

% Shape to a matrix with TEs second and all the rest first
vectorReorderDecoded = reshape(vectorReorderDecoded, size(decodedPWI4D, 4)/nInnerRepetitions, nInnerRepetitions);

% Flip the two dimensions and make a row vector again
vectorReorderDecoded = reshape(vectorReorderDecoded', 1, size(decodedPWI4D, 4));

% Reorder the data
decodedPWI4D = decodedPWI4D(:,:,:,vectorReorderDecoded);



%% 9. Calculate the correct parameters of the decoded data
% Calculate the correct parameters
indexDecoding = ((nInnerRepetitions+1):nVolumesPerInnerRepetition)'*ones(1,nAcquisitionRepetitions) + ones(nVolumesPerInnerRepetition-nInnerRepetitions,1)*(0:nAcquisitionRepetitions-1)*nVolumesPerInnerRepetition;
indexDecoding = indexDecoding(:)';
xQ.EchoTime_PWI4D = xQ.EchoTime(indexDecoding);
xQ.InitialPLD_PWI4D = xQ.Initial_PLD(indexDecoding);
xQ.LabelingDuration_PWI4D = xQ.LabelingDuration(indexDecoding);
    

end