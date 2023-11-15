function [x] = xASL_vis_AddIM2QC(x,parms)
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
% Copyright 2015-2023 ExploreASL

    %% Admin
	if nargin < 2
		error('Two parameters required');
	end

	parms = xASL_HandleInputPars(parms, 'bCrop', true); % crop by default
    parms = xASL_HandleInputPars(parms, 'FileName', 'n/a');
    [parms, DidntContain] = xASL_HandleInputPars(parms,'IM',[]);
    
    if  DidntContain
        warning('Parms.IM (input image) missing, aborting');
        return;
    else
        IM  = parms.IM;
    end
    
    [parms, DidntContain] = xASL_HandleInputPars(parms,'ModuleName',[]);    
    
    if  DidntContain
        warning('Parms.ModuleName missing, aborting');
        return;
    end    
    
    if  parms.bCrop
        X = size(IM,1); Y = size(IM,2);
        IM = squeeze(IM(ceil(0.33*X)+2:floor(0.67*X)-1,ceil(Y/4+1):floor(Y/2),:)); % slice 6
    end

    if  xASL_stat_SumNan(IM(:))==0 % if the image was empty
        return; % exit function
    end    

    %% Create the fields
    if ~isfield(x,'Output_im') || isempty(fields(x.Output_im))
        IndexIm = 1;
        x.Output_im = struct;
	end
    if ~isfield(x.Output_im, parms.ModuleName)
        IndexIm = 1;
	elseif strcmpi(parms.ModuleName, 'structural')
		if ~iscell(x.Output_im.(parms.ModuleName))
			x.Output_im = rmfield(x.Output_im, parms.ModuleName);
			IndexIm = 1;
		else
			IndexIm = length(x.Output_im.(parms.ModuleName)) + 1;
		end
	else
		if ~isfield(x.Output_im.(parms.ModuleName), x.SESSION)
			IndexIm = 1;
		elseif ~iscell(x.Output_im.(parms.ModuleName).(x.SESSION))
			x.Output_im.(parms.ModuleName) = rmfield(x.Output_im.(parms.ModuleName), x.SESSION);
		else
			IndexIm = length(x.Output_im.(parms.ModuleName).(x.SESSION)) + 1;
		end
    end
        
    %% Add the image to the field
	if strcmpi(parms.ModuleName, 'structural')
		x.Output_im.(parms.ModuleName){IndexIm} = IM;
	else
		x.Output_im.(parms.ModuleName).(x.SESSION){IndexIm} = IM;
	end

end


%% ========================================================================================
%% ========================================================================================
function [StructIn DidntContain] = xASL_HandleInputPars(StructIn,FieldName,DefaultV)
%xASL_HandleInputPars Summary of this function goes here
%   Detailed explanation goes here

DidntContain        = false;
if     ~isfield(StructIn,FieldName)
        DidntContain    = true;
elseif  isempty(StructIn.(FieldName))
        DidntContain    = true;
end

if  DidntContain
    StructIn        = setfield(StructIn,FieldName,DefaultV); % create field with default value
end


end