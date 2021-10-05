function [x] = xASL_imp_DetermineStructureFromRawdata(x)
%xASL_imp_DetermineStructureFromRawdata Determine structure from rawdata
%
% FORMAT: [x] = xASL_imp_DetermineStructureFromRawdata(x)
%
% INPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Determine structure from rawdata.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% SUBJECTS
    % x.SUBJECTS = xASL_adm_GetFileList(x.modules.import.imPar.TempRoot,[],false,[],true);
    
    %if isempty(x.SUBJECTS)
    %    warning('Unable to find subjects in temp directory...');
    %end
    
    error('Implement this...');
    


end


