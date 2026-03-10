function [x] = xASL_init_LoadDataPar(x)
%xASL_init_LoadDataPar Load data parameter file
%
% FORMAT: [x] = xASL_init_LoadDataParameterFile(x)
%
% INPUT:
%   x       - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:
%   x       - ExploreASL x structure (STRUCT)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Load data parameter file (which contains settings for running the pipeline, 
%              so should actually be called something like xASLProcessingSettings.json
%
% EXAMPLE:     This is part of the initialization workflow. Check out the usage there.
%
% This function performs the following parts:
% 1. Generate warning for incompatibility with old dataPar.m
% 2. Choose the dataPar location
% 3. Load pre-existing dataPar
% 4. Populate dataPar with missing parameters
% 5. Write final dataPar
% 6. Load dataPar
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:  n/a
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

    
    

    %% 1. Generate warning for incompatibility with old dataPar.m
    [~, ~, Dext] = fileparts(x.dir.dataPar);
    if strcmp(Dext,'.m')
        warning('No .m file backwards compatibility starting v1.10.0...');
    elseif strcmp(Dext,'.mat')
        warning('No .mat file backwards compatibility starting v1.10.0...');
    end


    %% 2. Choose the dataPar location
    bUseRoot = false;
    bUseRawdata = false;
    bUseDerivatives = false;

    % Check dataPar inside the root folder /
    listRoot = xASL_adm_GetFileList(x.dir.DatasetRoot, '(?i)(^dataPar.*\.json$)', 'FPList', [], 0);
    if ~isempty(listRoot)
        bUseRoot = true;
    end   

    % Check dataPar inside the rawdata folder /rawdata
    listRawdata = xASL_adm_GetFileList(x.dir.RawData, '(?i)(^dataPar.*\.json$)', 'FPList', [], 0);
    if ~isempty(listRawdata)
        warning([xASL_num2str(length(listRawdata)) ' dataPar.json (or similar) file(s) found in ' x.dir.RawData ', will try this but this is not the appropriate location']);
        bUseRawdata = true;
    end

    % Check dataPar inside the processing folder /derivatives/ExploreASL
    fListLegacy = xASL_adm_GetFileList(x.dir.xASLDerivatives, '(?i)(^dataPar.*\.json$)', 'FPList', [], 0);
    if ~isempty(fListLegacy)
        bUseDerivatives = true;
    end

    % Choose folder & filelist
    if bUseDerivatives && bUseRoot
            % We prioritize existing dataPar.json in the derivatives folder
        fprintf('%s\n', 'dataPar.json (or similar) file(s) both in /derivatives/ExploreASL and / root folder, ignoring root folder');
        bUseRoot = false;
    end
    if bUseDerivatives
        fFolder = x.dir.xASLDerivatives;
        listDatapar = fListLegacy;
        fprintf('%s\n', ['dataPar.json found in ' fFolder]);
    elseif bUseRoot
        fFolder = x.dir.DatasetRoot;
        listDatapar = listRoot;
        fprintf('%s\n', ['dataPar.json found in ' fFolder]);
    elseif bUseRawdata
        fFolder = x.dir.RawData;
        listDatapar = listRawdata;
        fprintf('%s\n', ['dataPar.json found in ' fFolder]);
    else
        listDatapar = {};
    end


    %% 3. Load pre-existing dataPar
    if length(listDatapar)>1
        fprintf('Warning: multiple dataPar.json files found, using the first\n');
        fprintf('%s\n', [': ' listDatapar{1}]);
        dataPar.x = xASL_io_ReadDataPar(listDatapar{1});
    elseif isempty(listDatapar)
        % Create default if no dataPar was provided
        fprintf('No dataPar.json provided, will create default version\n');
        dataPar = struct();
    else
        dataPar.x = xASL_io_ReadDataPar(listDatapar{1});
    end


    %% 4. Populate dataPar with missing parameters
    % Fill in important information in the dataPar if missing
    if ~isfield(dataPar, 'x')
        % Add x field
        dataPar.x = struct;
    end

    % Check for settings fields
    if ~isfield(dataPar.x, 'settings')
        dataPar.x.settings = struct;
    end
    % Check for quality field
    if ~isfield(dataPar.x.settings, 'Quality')
        dataPar.x.settings.Quality = true;
    end
    % Check for DELETETEMP field
    if ~isfield(dataPar.x.settings, 'DELETETEMP')
        dataPar.x.settings.DELETETEMP = true;
	end

	% Check for Atlases and TissueMasking parameters
	dataPar.x = xASL_init_LoadDataPar_PrepareAtlas4ROI(dataPar.x);


    %% 5. Write final dataPar
    x.dir.dataPar = fullfile(x.dir.xASLDerivatives, 'dataPar.json');
    xASL_io_WriteJson(x.dir.dataPar, dataPar);

    
    %% 6. Load dataPar
    x = xASL_adm_MergeStructs(dataPar.x, x);
    

end


%% -----------------------------------------------------------------------------
%% -----------------------------------------------------------------------------
function [x] = xASL_init_LoadDataPar_PrepareAtlas4ROI(x)
%Check for x.S.Atlases & x.S.TissueMasking validity and print instructions if incorrect values are provided

bAtlasTissueMatch = true; % x.S.Atlases & x.S.TissueMasking should match before we can continue (needed for the population module)

if isfield(x, 'S') % Check upfront because we don't want to create a potentially empty S field without subfields
	% First ensure that these fields are cells
    if isfield(x.S, 'Atlases') && ~iscell(x.S.Atlases)
        x.S.Atlases = {x.S.Atlases};
    end
    if isfield(x.S, 'TissueMasking') && ~iscell(x.S.TissueMasking)
        x.S.TissueMasking = {x.S.TissueMasking};
    end
    
    if ~isfield(x.S,'Atlases') && ~isfield(x.S, 'TissueMasking')
		% Default atlases/ROIs & tissue masks are provided in the population module
	elseif ~isfield(x.S, 'Atlases') && isfield(x.S, 'TissueMasking')
		% Missing Atlases, but provided TissueMasking - cannot continue
		warning('Custom tissue-types (x.S.TissueMasking) specified without ROI atlas-selection (x.S.Atlases). Atlases need to be provided. See instructions below:');
		bAtlasTissueMatch = false;
	elseif isfield(x.S, 'Atlases') && isfield(x.S, 'TissueMasking') && length(x.S.Atlases)~=length(x.S.TissueMasking)
		% Non matching lengths, cannot continue
		warning('The number of ROI atlases x.S.Atlases as subject-wise tissue-types x.S.TissueMasking provided does not match:');
		fprintf('%s\n', ['S:{Atlases:["' strjoin(x.S.Atlases, '", "') '"]}']);
		fprintf('%s\n', ['S:{TissueMasking:["' strjoin(x.S.TissueMasking, '", "') '"]}']);
		bAtlasTissueMatch = false;
	elseif 	isfield(x.S, 'Atlases') && ~isfield(x.S, 'TissueMasking')
		% TissueMasking not provided, so it has to be extracted from Atlases as previously.
		warning('ROIs provided in S.Atlases without the tissue-types for these ROIs in S.TissueMasking.');
		bAtlasTissueMatch = false;
		textRecommendation = 'S:{TissueMasking:["';
		% Fill in TissueMasking based on the atlas names and default to GM. The user will see a warning and automatic tissue masks to verify
		for iAtlas = 1:numel(x.S.Atlases)
			if iAtlas > 1
				textRecommendation = [textRecommendation '", "'];
			end
			if ~isempty(regexpi(x.S.Atlases{iAtlas}, 'WM')) || ~isempty(regexpi(x.S.Atlases{iAtlas}, 'whitematter'))
				textRecommendation = [textRecommendation 'WM'];
			elseif ~isempty(regexpi(x.S.Atlases{iAtlas}, 'WB')) || ~isempty(regexpi(x.S.Atlases{iAtlas}, 'wholebrain'))
				% In case the atlas name contains WB, we inform the user that the probable tissue mix the user wants is GM+WM. This is a recommendation only and the user needs to specify in the configuration if GM+WM or GM+WM+CSF is wanted
				textRecommendation = [textRecommendation 'GM+WM'];
			else
				textRecommendation = [textRecommendation 'GM'];
			end
		end
		textRecommendation = [textRecommendation '"]}'];
		fprintf('Recommended addition to the dataPar.json to match provided Atlases: %s\n\n', textRecommendation);
	end

	if isfield(x.S, 'TissueThreshold')
		if ~isfield(x.S,'Atlases') && ~isfield(x.S, 'TissueMasking') && (length(x.S.TissueThreshold) > 1)
			% If TissueThreshold is provided without the other two parameters, it has to have length 1
			warning('TissueThreshold has to have one value only when TissueMasking and Atlases are not provided');
		elseif (isfield(x.S, 'Atlases') && length(x.S.Atlases)~=length(x.S.TissueThreshold)) || (isfield(x.S, 'TissueMasking') && length(x.S.TissueMasking)~=length(x.S.TissueThreshold))
			warning('Length of Atlases and TissueMasking has to match the length of TissueThreshold');
		end
	end

end
if ~bAtlasTissueMatch
    fprintf('%s\n', 'When ROI atlases are provided in S.Atlases, their tissue types,');
    fprintf('%s\n', 'you need to be provided as well in dataPar.json as "S":{"TissueMasking"}, with either option "GM", "WM", "GM+WM".');
	warning('x.S.TissueMasking input did not match x.S.Atlases');
	fprintf('%s\n', 'If you want GM CBF from an atlas ROI, add "GM". If you want WM CBF from an atlas ROI, add "WM"');
	fprintf('%s\n', 'If you want both GM CBF and WM CBF from an atlas ROI, add the atlas ROI twice and add both "GM" and "WM"');
	fprintf('%s\n', 'E.g., if you want both GM CBF and WM CBF from the Harvard-Oxford atlas, you would provide the following in dataPar.json:');
	fprintf('%s\n', '"S":{"Atlases":["HOConn", "HOConn"]}');
	fprintf('%s\n', '"S":{"TissueMasking":["GM", "WM"]}');

	% No match means that we have to end it
    error('Not the same number of ROI atlases as subject-wise tissue-types, skipping');
end


end
