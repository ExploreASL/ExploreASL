function jsonOut = xASL_bids_BIDSifyASLNII(jsonIn, bidsPar, pathIn, pathOutPrefix)
%xASL_bids_BIDSifyASLNII Takes the NII file and saves into the correct NII format
%
% FORMAT: jsonOut = xASL_bids_BIDSifyASLNII(jsonIn, bidsPar, pathIn, pathOutPrefix)
%
% INPUT:
%   jsonIn        - JSON with the input fields in BIDS format (REQUIRED)
%   bidsPar       - Structure with the BIDS configuration information (REQUIRED)
%   pathIn        - Path to the source ASL NIfTI file (REQUIRED)
%   pathOutPrefix - Path to the output - a prefix only that needs to be extended by ASL.nii or aslContext.tsv (REQUIRED)
%
% OUTPUT: 
%   jsonOut       - Output JSON, the scalings and ASLContext is removed
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:     It modifies the NIfTI file to take into account several BIDS
%                  specifics. Specifically, it applies the previously calculated
%                  scalings, and  it saves the ASLcontext.tsv file, 
%
% EXAMPLE: n/a
%
% __________________________________
% Copyright 2015-2021 ExploreASL

% Remove the AslContext field and save it as a separate file
filenameTSV = [pathOutPrefix '_' bidsPar.stringAslContext '.tsv'];
[pathTSV,~,~] = fileparts(filenameTSV);

xASL_adm_CreateDir(pathTSV);

% Write volume types (asl context) to file
fContext = fopen(filenameTSV,'w+');
fwrite(fContext,sprintf('volume_type\n'));
fwrite(fContext,jsonIn.ASLContext);
fclose(fContext);

% Remove context field
jsonOut = rmfield(jsonIn,'ASLContext');

% Validate ASL NIFTI output path
[~,outputFileASL,outputExtensionASL] = xASL_fileparts([pathOutPrefix '_asl.nii.gz']);
outputFilenameASL = [outputFileASL outputExtensionASL];
xASL_bids_ValidateNiftiName(outputFilenameASL,'asl');
	
% Save the JSON and NII to final location for ASL
headerASL = xASL_io_ReadNifti(pathIn);
if (jsonOut.scaleFactor && jsonOut.scaleFactor~=1) || length(headerASL.dat.dim) < 4 || headerASL.dat.dim(4) == 1
    imNii = xASL_io_Nifti2Im(pathIn);
    
    if (jsonOut.scaleFactor && jsonOut.scaleFactor~=1)
        imNii = imNii .* jsonOut.scaleFactor;
    end
    
    % Scaling changed, so we have to save again OR
    % The fourth dimension is 1, so we have to write the file again, to make sure the
    xASL_io_SaveNifti(pathIn,[pathOutPrefix '_asl.nii.gz'],imNii,[],1,[]);
    
    % Delete original Nifti if the in and out file have different paths
    if ~strcmp(pathIn,[pathOutPrefix '_asl.nii.gz']) && ~strcmp(pathIn,[pathOutPrefix '_asl.nii'])
        xASL_delete(pathIn);
    end
else
	% Move the ASL
	xASL_Move(pathIn,[pathOutPrefix '_asl.nii.gz'],1);
end

jsonOut = rmfield(jsonOut,'scaleFactor');

end


