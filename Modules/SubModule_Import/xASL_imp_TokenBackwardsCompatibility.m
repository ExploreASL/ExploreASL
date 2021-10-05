function imPar = xASL_imp_TokenBackwardsCompatibility(imPar)
%xASL_imp_TokenBackwardsCompatibility Fix the token backwards compatibility
%
% FORMAT: imPar = xASL_imp_TokenBackwardsCompatibility(imPar)
%
% INPUT:
%   imPar - JSON file with structure with import parameters (REQUIRED, STRUCT)
%
% OUTPUT:
%   imPar - JSON file with structure with import parameters
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Fix the token backwards compatibility.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
%
% __________________________________
% Copyright 2015-2021 ExploreASL



    % Check token ordering
    if length(imPar.tokenOrdering)==3
        % Backwards compatibility visits
        if size(imPar.tokenOrdering,1) > 1
            % Vertical vector
            imPar.tokenOrdering = imPar.tokenOrdering';
        end
        imPar.tokenOrdering = [imPar.tokenOrdering(1) 0 imPar.tokenOrdering(2:3)]; % insert Visits(none)
    elseif length(imPar.tokenOrdering)==2 
        % Backwards compatibility Visits & sessions
        imPar.tokenOrdering = [imPar.tokenOrdering(1) 0 0 imPar.tokenOrdering(2)]; % insert sessions & visits(none)
    end


    % Change dcmnii_version for PARREC if needed
    if ~isempty(strfind(char(imPar.folderHierarchy(end)),'PAR'))
        imPar.dcm2nii_version = '20101105';
    end


end


