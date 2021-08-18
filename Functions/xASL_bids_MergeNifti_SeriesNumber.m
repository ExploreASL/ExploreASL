function pathOut = xASL_bids_MergeNifti_SeriesNumber(NiftiPaths, niiTable)
%xASL_bids_MergeNifti_SeriesNumber Take a list of NIfTI files and
%concatenates 3D/4D files into a 4D sequence if possible according to SeriesNumber for niiTable
%
% FORMAT: pathOut = xASL_bids_MergeNifti_SeriesNumber(NiftiPaths, niiTable)
% 
% INPUT:
%   NiftiPaths - cell containing list of strings with full paths of the files (REQUIRED)
%   niiTable   - cell containing a table Filename, InstanceNumber, SeriesNumber, FileType, FilePath (REQUIRED)
% OUTPUT:
%   pathOut    - output paths
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Take a list of NIfTI files and concatenates 3D/4D files into a 4D sequence if possible according to SeriesNumber for niiTable
%
% EXAMPLE:     pathOut = xASL_bids_MergeNifti_SiemensASLFiles(NiftiPaths, niiTable);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2021 ExploreASL

% 0. Admin
% Nothing merged
pathOut = '';

% Checks if the niiTable exists
if nargin < 2 || isempty(niiTable)
	return;
end

% The SeriesNumber column exist and is numeric
if size(niiTable,2) < 5 || ~isnumeric(niiTable{1,5})
	return;
end

% All cells are filled 
if sum(cellfun('isempty',niiTable(:,5))) > 0
	return;
end

% Initialize the vector
vectorSeriesNumber = zeros(size(niiTable,1),1);
for iField = 1:size(niiTable,1)
	vectorSeriesNumber(iField,1) = niiTable{iField,5}(1);
end

if sum(isnan(vectorSeriesNumber)) > 0
	return;
end

% Sort and extract unique numbers
[~,indexSort,~] = unique(vectorSeriesNumber);

% In case all numbers are unique, then sort the files accordingly
if length(indexSort) == length(vectorSeriesNumber)
	pathOut = xASL_bids_MergeNifti_Merge(NiftiPaths,indexSort,'ASL4D',0);
end

end
