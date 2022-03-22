function [x] = xASL_init_SubStructs(x)
%xASL_init_SubStructs Initialize ExploreASL x structure substructs.
%
% FORMAT: [x] = xASL_init_SubStructs(x)
%
% INPUT:
%   x       - ExploreASL x struct (STRUCT, REQUIRED)
%
% OUTPUT:
%   x       - ExploreASL x struct
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Initialize the ExploreASL x structure substructs/fields.
%               Only fields which do not exist so far are added.
%               This script is supposed to help with the overall modularity of ExploreASL.
%               This script is identical to the function ExploreASL_Initialize_SubStructs within
%               ExploreASL_Initialize. We can not call this script from ExploreASL_Initialize,
%               since the paths are not initialized at that part of the script yet.
%
% EXAMPLE:      [x] = xASL_init_SubStructs(x);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright (c) 2015-2022 ExploreASL

    if nargin < 1 || isempty(x)
		x = struct;
	end
    
    % Statistics, directories, paths, and sequence related fields
    if ~isfield(x,'S'),                     x.S = struct;                   end
    if ~isfield(x,'D'),                     x.D = struct;                   end
    if ~isfield(x,'P'),                     x.P = struct;                   end
    if ~isfield(x,'Q'),                     x.Q = struct;                   end
    
    % Module subfields (import, structural, asl, & population)
    if ~isfield(x,'modules'),               x.modules = struct;             end
    if ~isfield(x.modules,'import'),        x.modules.import = struct;      end
    if ~isfield(x.modules,'structural'),    x.modules.structural = struct;  end
    if ~isfield(x.modules,'asl'),           x.modules.asl = struct;         end
    if ~isfield(x.modules,'population'),    x.modules.population = struct;  end
    
    % Dataset related fields, workflow settings, toolbox/external (SPM, CAT, FSL, etc.) fields, general directories
    if ~isfield(x,'dataset'),               x.dataset = struct;             end
    if ~isfield(x,'settings'),              x.settings = struct;            end
    if ~isfield(x,'external'),              x.external = struct;            end
    if ~isfield(x,'dir'),                   x.dir = struct;                 end     
    if ~isfield(x,'opts'),                  x.opts = struct;                end     

end

