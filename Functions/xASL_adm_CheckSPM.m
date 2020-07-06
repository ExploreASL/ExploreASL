function [spm_path, spm_version] = xASL_adm_CheckSPM(modality, proposed_spm_path, check_mode)
% Checks if the spm function exists and if the reported version matches our development version.
%
% FORMAT: [spm_path, spm_version] = xASL_adm_CheckSPM([modality, proposed_spm_path, check_mode])
%         [spm_path]              = xASL_adm_CheckSPM(...)
%                                   xASL_adm_CheckSPM(...)
%
% INPUT:
%   modality          - Modality of SPM check. Can be 'FMRI' or 'FULL' (default = 'FMRI')
%   proposed_spm_path - Is the user specified path to the SPM toolbox (no default) 
%   check-mode        - Mode of SPM. Can be 'basic' or 'full' (default = 'basic')
% OUTPUT:
%   spm_path          - Path to the SPM folder
%   spm_version       - Version of SPM
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Checks if the spm function exists and if the reported version matches our development 
%              version ({{SPM8}} or {{SPM12}}). If the spm toolbox is not available yet, it will try the 
%              PROPOSED_SPM_PATH (if specified) or the user selected directory and add it to {{PATH}}.
%              The function will fail if {{SPM}} cannot be found or if detecting an unsupported version.
% EXAMPLE: [spm_path, spm_version] = xASL_adm_CheckSPM();
%          [spm_path] = xASL_adm_CheckSPM('FMRI');  
%          [spm_path, spm_version] = xASL_adm_CheckSPM('FMRI','/usr/local/spm12');
%          [spm_path, spm_version] = xASL_adm_CheckSPM('FMRI','/usr/local/spm12','basic');
%          [spm_path, spm_version] = xASL_adm_CheckSPM('FMRI','/usr/local/spm12','full');
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright © 2015-2020 ExploreASL

    % Check the input parameters
	if nargin<1 || isempty(modality)
		modality = 'FMRI';
	end

	if nargin<3 || isempty(check_mode)
		check_mode = 'basic';
	end
	
    % Check if SPM toolbox is already available
    p = which('spm');
	if isempty(p)
		if ~isdeployed % Adding in deployed versions does not make sense
			if nargin>=2 && isdir(proposed_spm_path)
				p = proposed_spm_path;
			else
				error('xASL_adm_CheckSPM:Can''t continue without SPM toolbox');
			end
			if ~exist(fullfile(p,'spm.m'),'file')
				error('xASL_adm_CheckSPM:Specified directory is not an SPM toolbox directory: [%s]',p);
			end
			addpath(p);
		end
	else
		% SPM already available, just strip /spm.m to return the path
		p = fileparts(p);
	end
	
    % Check SPM installation
    spm_check_installation(check_mode);

    % Check SPM version
    v = spm('ver');
	if ~strcmp(v,'SPM8') && ~strcmp(v,'SPM12')
		error('xASL_adm_CheckSPM:Compatible only with SPM8 and SPM12. Current SPM version [%s]',v);
	end
    
    % Configure SPM and batch system
    spm('defaults', modality);    
    spm_jobman('initcfg');
    
    % Return path and version if requested
    if nargout>=1
        spm_path = p;
    end
    if nargout>=2
        spm_version = v;
    end
end

