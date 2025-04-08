function [x] = xASL_vis_AddIM2QC(x, parms)
%xASL_vis_AddIM2QC Checks which images already are loaded, and  adds new image.
%
% FORMAT:       [x] = xASL_vis_AddIM2QC(x,parms);
% 
% INPUT:        
%       x        x-struct (REQUIRED)
%       parms    A structure containing parameters 
%                bCrop (DEFAULT true)
%                FileName (DEFAULT n/a)
%                IM (REQUIRED)
%                ModuleName - 'ASL' or 'Structural' (REQUIRED)
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Checks which images already are loaded, and  adds new image.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2024 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________


    %% Admin
	if nargin < 2
		error('Two parameters required');
	end

	parms = xASL_HandleInputPars(parms, 'bCrop', true); % crop by default
    parms = xASL_HandleInputPars(parms, 'FileName', 'n/a');
    
    if ~isfield(parms, 'IM') || isempty(parms.IM)
        warning('parms.IM (input image) missing, aborting');
        return;
    elseif xASL_stat_SumNan(parms.IM(:))==0 % if the image was empty
        return; % exit function
    else
        IM  = parms.IM;
    end

    if ~isfield(parms, 'ModuleName') || isempty(parms.ModuleName)
        warning('Parms.ModuleName missing, aborting');
        return;
    end    

    if ~isfield(parms, 'paths') || isempty(parms.paths)
        warning('parms.paths (input image filename) missing, aborting');
        return;
    end

    if  parms.bCrop
        X = size(IM,1); Y = size(IM,2);
        IM = squeeze(IM(ceil(0.33*X)+2:floor(0.67*X)-1,ceil(Y/4+1):floor(Y/2),:)); % slice 6
    end


    %% Create the field
    if ~isfield(x, 'Output_im') || isempty(x.Output_im)
        x.Output_im = struct;
	end


    %% Create fieldname from filename
    if ischar(parms.paths)
        pathsAre = {parms.paths};
    elseif iscell(parms.paths)
        pathsAre = parms.paths;
    else
        error('Unknown format of parms.paths');
    end

    if isfield(parms, 'preFix') && ~isempty(parms.preFix)
        pathName = parms.preFix;
        if ~strcmp(pathName(end), '_')
            pathName = [pathName '_'];
        end
    else
        pathName = [];
    end

    for iPath=1:length(pathsAre)
        [~, fFileName] = xASL_fileparts(pathsAre{iPath}); % filename
        fFileName = strrep(fFileName, x.SUBJECT, ''); % remove subjectname
        if isfield(x, 'SESSION') && ~isempty(x.SESSION)
            fFileName = strrep(fFileName, x.SESSION, ''); % remove sessionname
        end
        
        % Remove trailing underscore(s)
        [iStart, iEnd] = regexp(fFileName, '_*');
        
        % Take last index
        iStart = iStart(end);
        iEnd = iEnd(end);
        if iEnd==length(fFileName)
            fFileName = fFileName(1:iStart-1);
        end

        % Add overlays
        if iPath~=length(pathsAre)
            fFileName = [fFileName '_with_'];
        end
        pathName = [pathName fFileName];
    end
    

    %% Add the image to the field
	if strcmpi(parms.ModuleName, 'structural')
		%% ==============================================
		%% PM: BACKWARD COMPATIBILITY CODE, CAN BE PHASED OUT
		%% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
		% Until version 1.11.0 and 1.12.0_beta, images were added as cells. 
		% From 2.0.0, new images are added in the structure under a specific name.
		% When processing new subjects, the previous QC structure is loaded, so in case the old version 
		% is loaded, it has to be removed as it is not fully compatible
		if isfield(x.Output_im, parms.ModuleName) && iscell(x.Output_im.(parms.ModuleName))
			x.Output_im.(parms.ModuleName) = [];
		end
		%% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		%% PM: BACKWARD COMPATIBILITY CODE, CAN BE PHASED OUT
		%% ==============================================
		
		x.Output_im.(parms.ModuleName).(pathName) = IM;
	else
		%% ==============================================
		%% PM: BACKWARD COMPATIBILITY CODE, CAN BE PHASED OUT
		%% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
		if isfield(x.Output_im, parms.ModuleName) && isfield(x.Output_im.(parms.ModuleName), x.SESSION) && iscell(x.Output_im.(parms.ModuleName).(x.SESSION))
			% Fix compatibility issue with previously saved data
			x.Output_im.(parms.ModuleName).(x.SESSION) = [];
		end
		%% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		%% PM: BACKWARD COMPATIBILITY CODE, CAN BE PHASED OUT
		%% ==============================================
		x.Output_im.(parms.ModuleName).(x.SESSION).(pathName) = IM;
	end

end


%% ========================================================================================
%% ========================================================================================
function [StructIn, DidntContain] = xASL_HandleInputPars(StructIn, FieldName, DefaultValue)
%xASL_HandleInputPars Summary of this function goes here
%   Detailed explanation goes here

DidntContain = false;
if ~isfield(StructIn, FieldName) || isempty(StructIn.(FieldName))
    DidntContain = true;
end

if DidntContain
    StructIn = setfield(StructIn, FieldName, DefaultValue); % create field with default value
end


end