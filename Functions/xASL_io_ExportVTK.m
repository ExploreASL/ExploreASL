function xASL_io_ExportVTK(pathExploreASL,nifti,mask,exportPath)
%xASL_io_ExportVTK Export VTK image file.
%
% FORMAT: xASL_io_ExportVTK(pathExploreASL, nifti, [mask, exportPath])
%
% INPUT:
%   pathExploreASL - Path to ExploreASL (REQUIRED, CHAR ARRAY)
%   nifti          - Path to NIFTI image or image matrix (REQUIRED, CHAR ARRAY or 3D IMAGE)
%   mask           - Path to NIFTI mask (OPTIONAL, CHAR ARRAY, DEFAULT = [])
%   exportPath     - Path of the exported VTK file (OPTIONAL, DEFAULT = nifti path -> export.vtk)
%
%   OUTPUT:
%   n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      Export a VTK image file based on a 3D NIFTI or a 3D image matrix.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:          To export the "test" NIFTI to structured points in VTK format, you can run the following lines.
%                   [x] = ExploreASL_Initialize('',0);
%                   nifti = '.\test.nii';
%                   xASL_io_ExportVTK(x.MyPath,nifti);
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Check ExploreASL
    if nargin < 1
        warning('Could not find ExploreASL...');
        return 
    end

    %% Load image
    if nargin < 2
       warning('Neither input image path nor matrix given...');
       return 
    end
    
    if ischar(nifti)
        % Import NIFTI
        image = xASL_io_Nifti2Im(nifti);
    elseif isnumeric(nifti)
        % Import image
        image = nifti;
    end
    
    % Check dimensions
    if numel(size(image))~=3
        error('Input image dimensions are incorrect...');
    end

    % Apply mask
    if nargin > 2 && exist('mask','var')
        if ~isempty(mask)
            % Get mask
            imageMask = xASL_io_Nifti2Im(mask);
            % Resample mask to image
            [maskResampled] = xASL_im_ResampleLinearFair(imageMask,size(image));
            % Apply mask
            image = image.*maskResampled;
        end
    end
    
    % Check export path
    if nargin < 4 && ~exist('exportPath','var')
        if ischar(nifti)
            % Define export file
            folder = fileparts(nifti);
            exportPath = fullfile(folder,'export.vtk');
        else
            % Neither export path nor nifti path given
            exportPath = fullfile(pwd,'export.vtk');
        end
    end
    
    % Select export method
    fprintf('Export vtk as structured points...\n');
    image = squeeze(image);
    % Export to vtk (MIT-License) % https://www.mathworks.com/matlabcentral/fileexchange/47814-vtkwrite-exports-various-2d-3d-data-to-paraview-in-vtk-file-format
    readAndPrintTextFile(fullfile(pathExploreASL,'External','vtkwrite','MIT_License.txt'));
    vtkwrite(exportPath, 'structured_points', 'image', image)

end


%% Read text file into cell array and print it
function readAndPrintTextFile(filePath)

    % Display which license is being printed
    fprintf('================================================================================\n')
    folderStructure=regexp(filePath,filesep,'split');
    fprintf('%s\n',folderStructure{end-1},1);

    % Read file
    fileStr = fileread(filePath);
    fileLines = regexp(fileStr, '\r\n|\r|\n', 'split');
    fileLines = fileLines';
    
    % Print file
    for thisLine=1:numel(fileLines)
        fprintf('%s\n',fileLines{thisLine,1});
    end
    fprintf('================================================================================\n')

end


