function xASL_wrp_CreatePopulationTemplates(x, bSaveUnmasked, bCompute4Sets, SpecificScantype, bSkipWhenMissingScans, bRemoveOutliers, FunctionsAre, bUpdateMetadata, SmoothingFWHM, bSaveMasks4QC)
%xASL_wrp_CreatePopulationTemplates ExploreASL Population module wrapper,
%creates population parametric images for each ScanType
%
% FORMAT: xASL_wrp_CreatePopulationTemplates(x[, bSaveUnmasked, Compute4Sets, SpecificScantype, bSkipWhenMissingScans, bRemoveOutliers, FunctionsAre])
%
% INPUT:
%   x            - structure containing fields with all information required to run the population module (REQUIRED)
%   x.bForceTemplates - boolean to force creation templates, even with a very
%                       low number of subjects
%   bSaveUnmasked - allows saving the same images without masking (OPTIONAL, DEFAULT=true)
%   bCompute4Sets - creates the same parametric images for subsets of the
%                  data (e.g. cohorts, sites, etc)
%                  can be a Boolean for computing for all subsets (1) or
%                  none (0), or a cell structure with the names of the
%                  subsets to compute for. (OPTIONAL, DEFAULT=0)
%   SpecificScantype - allows providing ScanTypes to create templates for:
%                      should be 3 cells:
%                      SpecificScantype{1} = PreFixList (prefix string of the NIfTI, associated with the ScanType, e.g. 'qCBF', 'PV_PGM', etc)
%                      SpecificScantype{2} = TemplateNameList (output name to give to the NIfTI file of the template)
%                      SpecificScantype{3} = SessionsExist (whether sessions exist, if true it uses ASL_1 ASL_2 ASL_n)
%                      
%                      (OPTIONAL, DEFAULT = use predefined ScanTypes)
%   bSkipWhenMissingScans - This parameter allows to choose if we want to
%                       skip creating templates if more than 10% of the datasets are
%                       present (1), or also create templates when subjects
%                       are missing (0). 
%                       (OPTIONAL, DEFAULT=0)
%   bRemoveOutliers   - This parameter makes robust statistics, by removing
%                       outliers (OPTIONAL, DEFAULT=false)
%   FunctionsAre      - two cells, with functions, and with names
%                       (OPTIONAL, DEFAULT is mean and SD). Same setup as
%                       "SpecificScantype" above). First should be a
%                       functionhandle, e.g. @xASL_stat_MeanNan (not a
%                       string)
%   bUpdateMetadata   - boolean specifying if we reload the metadata, e.g.
%                       for potentially other defined cohorts etc in the
%                       participants.tsv. This can take some time though.  Only relevant when computing multiple sets. 
%                       (OPTIONAL, DEFAULT = false);
%   SmoothingFWHM     - Full-Width-Half-Maximum in [X Y Z] voxels for smoothing of the output image
%                       (OPTIONAL, DEFAULT = [0 0 0] (i.e. no smoothing)
%   bSaveMasks4QC     - boolean specifying if we wish to save the final
%                       images used for computations. Only recommended for
%                       masks, to avoid space redundancy, which is also why
%                       this is stored in uint8 and zipped. (OPTIONAL,
%                       DEFAULT = 0)
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function creates simple parametric images, a.k.a. templates, for
% different image/scan types, on population level, as well as for different
% sets (e.g. sites/scanners/cohorts, etc) if specified. By default these
% images are masked, and transformed into a single column, for quick
% computations with low memory usage. The default parametric images that
% are created are the mean, between-subject SD, and the maximal intensity
% projection (MIP). The latter can e.g. identify intra-vascular signal that
% is similar between different subjects. Other parametric maps can be
% decommented (now commented out for speed).
%
% Any new addition to participants.tsv will be recognized and loaded, for
% the generation of new parametric maps for groups specifically
% (needs to be set in input argument bCompute4Sets)
%
% If a set only includes a combination of the following SetOptions:
% left, right, l, r, n/a, NaN (irrespective of capitals)
% each image with option right/r, will be flipped in the left-right
% direction, and left/right will not be treated as separate groups.
% This function performs the following steps:
% 
% 1. Define images/scantypes (if they are not defined by input argument SpecificScantype)
% 2. Iterate over scan types & sessions
% 3. Check availability images
% 4. Load images
% 5. Remove outliers
% 6. Compute templates for all subjects together (only for bilateral images)
% 7. Compute templates for individual sets
%
% EXAMPLE: xASL_wrp_CreatePopulationTemplates(x);
% EXAMPLE for specific scantypes:
%          xASL_wrp_CreatePopulationTemplates(x, [], [], {{'qCBF'} {'CBF'} 1});
% EXAMPLE for specific scantypes and specific function:
%          xASL_wrp_CreatePopulationTemplates(x, 0, 1, {'qCBF' 'CBF' 1}, 1, 0, {{@xASL_stat_MeanNan} {'mean'}});
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2023 ExploreASL


% ----------------------------------------------------------------------------------------------------
%%  0. Admin

if nargin<2 || isempty(bSaveUnmasked)
    bSaveUnmasked = true;
end

if nargin<3 || isempty(bCompute4Sets)
    bCompute4Sets = 0;
elseif iscell(bCompute4Sets)
    Sets2Check = bCompute4Sets;
    bCompute4Sets = 1;
elseif bCompute4Sets==1
    Sets2Check = [];
elseif bCompute4Sets==0
    Sets2Check = [];
    % we also set this when computing 4 Sets is disabled, just to be sure
else
    error('Invalid bComputeSets option, skipping');
end

if nargin<5 || isempty(bSkipWhenMissingScans)
	bSkipWhenMissingScans = false;
end

if nargin<6 || isempty(bRemoveOutliers)
    bRemoveOutliers = false;
end

if nargin<7 || isempty(FunctionsAre)
    FunctionsAre{1} = {@xASL_stat_MeanNan,@xASL_stat_StdNan};
    FunctionsAre{2} = {'mean' 'sd'};
else
    if ~iscell(FunctionsAre{1})
        FunctionsAre{1}{1} = FunctionsAre{1};
    end
    if ~iscell(FunctionsAre{2})
        FunctionsAre{2}{1} = FunctionsAre{2};
    end    
end 

if nargin<8 || isempty(bUpdateMetadata)
    bUpdateMetadata = false;
end

if nargin<9 || isempty(SmoothingFWHM)
    SmoothingFWHM = [0 0 0];
elseif length(SmoothingFWHM)~=3 || ~isnumeric(SmoothingFWHM)
    error('Incorrect size of SmoothingFWHM kernel input, should have three numerical values');
elseif any(SmoothingFWHM<0)
    error('Negative values in smoothing kernel are invalid');
else
    SmoothingFWHM = double(SmoothingFWHM);
end

if nargin<10 || isempty(bSaveMasks4QC)
    x.S.bSaveMasks4QC = 0;
elseif bSaveMasks4QC==1
    x.S.bSaveMasks4QC = 1;
elseif bSaveMasks4QC==0
    x.S.bSaveMasks4QC = 0;
else
    error('Illegal value for bSaveMasks4QC');
end

if ~isfield(x,'GradualSkull')
    x.GradualSkull = xASL_io_Nifti2Im(fullfile(x.D.MapsSPMmodifiedDir, 'rbrainmask.nii'));
end

if ~isfield(x, 'bForceTemplates') || isempty(x.bForceTemplates)
    x.bForceTemplates = false;
elseif x.bForceTemplates==1
    bSkipWhenMissingScans = false;
end

if ~x.bForceTemplates && (x.dataset.nSubjects * x.dataset.nSessions < 6)
    % With too small datasets, created templated won't be reliable
    fprintf('\n\n%s\n\n', ['Only n=' xASL_num2str(x.dataset.nSubjects * x.dataset.nSessions) ' subject*runs, group templates may not be useful']);
end

Size1 = sum(x.S.masks.WBmask(:));

x.D.TemplatesStudyDir = fullfile(x.D.PopDir, 'Templates');
xASL_adm_CreateDir(x.D.TemplatesStudyDir);

if bCompute4Sets
    % Reload set parameters to be sure
    if bUpdateMetadata
        x = xASL_init_LoadMetadata(x); % Add statistical variables, if there are new ones
    end
    
    if isempty(Sets2Check)
        % Define sets to create comparative stats for
        for iSet=1:size(x.S.SetsID,2)
            TempUnique = unique(x.S.SetsID(:,iSet));
            nUniqueValues = length(TempUnique(~isnan(TempUnique)));
            if  nUniqueValues>1 && nUniqueValues<6
                Sets2Check(end+1,1) = iSet;
            end
        end
    else
        % convert names of sets to indices
        FoundSets = 0;
        for iSetCheck=1:length(Sets2Check)
            TempCheck = find(cellfun(@(y) strcmp(Sets2Check{iSetCheck}, y), x.S.SetsName));
            if isempty(TempCheck)
                warning(['Couldnt find set: ' Sets2Check{iSetCheck}]);
            else
                TempCheck2(iSetCheck) = TempCheck;
                FoundSets = FoundSets+1;
            end
        end
        Sets2Check = TempCheck2;
        if FoundSets==0
            error('No sets found, skipping. Please check participants.tsv vs bCompute4Sets');
        end
    end
end


% ----------------------------------------------------------------------------------------------------
%%  1. Define images/scantypes (if they are not defined by input argument SpecificScantype)
UsePredefined = true;
if nargin>3 && ~isempty(SpecificScantype)
    if ~iscell(SpecificScantype)
        warning('SpecificScantype should be a cell, running all predefined ScanTypes instead');
    elseif length(SpecificScantype)~=3
        warning('SpecificScantype should have 3 cells, running all predefined ScanTypes instead');
    else
        UsePredefined = false;
        if iscell(SpecificScantype{1})
            PreFixList = SpecificScantype{1};
        else
            PreFixList{1} = SpecificScantype{1};
        end
        if iscell(SpecificScantype{2})
            TemplateNameList = SpecificScantype{2};
        else
            TemplateNameList{1} = SpecificScantype{2};
        end
        
        SessionsExist = SpecificScantype{3};
    end
end
       
if UsePredefined
    % Structural images
    PreFixList          = {['r' x.P.STRUCT]};                   TemplateNameList           = {'T1'};        SessionsExist        =  0;
    PreFixList{end+1}   = ['mrc1' x.P.STRUCT];                  TemplateNameList{end+1}    = 'mrc1T1';      SessionsExist(end+1) =  0;
    PreFixList{end+1}   = ['mrc2' x.P.STRUCT];                  TemplateNameList{end+1}    = 'mrc2T1';      SessionsExist(end+1) =  0;
    PreFixList{end+1}   = ['rc1' x.P.STRUCT];                   TemplateNameList{end+1}    = 'pGM';         SessionsExist(end+1) =  0;
    PreFixList{end+1}   = ['rc2' x.P.STRUCT];                   TemplateNameList{end+1}    = 'pWM';         SessionsExist(end+1) =  0;
    PreFixList{end+1}   = ['rc3' x.P.STRUCT];                   TemplateNameList{end+1}    = 'pCSF';        SessionsExist(end+1) =  0;
    PreFixList{end+1}   = ['r' x.P.FLAIR];                      TemplateNameList{end+1}    = 'FLAIR';       SessionsExist(end+1) =  0;
    PreFixList{end+1}   = ['r' x.P.WMH_SEGM];                   TemplateNameList{end+1}    = 'WMH_SEGM';    SessionsExist(end+1) =  0;
	PreFixList{end+1}   = ['r' x.P.T1c];                        TemplateNameList{end+1}    = 'T1c';         SessionsExist(end+1) =  0;
	PreFixList{end+1}   = ['r' x.P.T2];                         TemplateNameList{end+1}    = 'T2';          SessionsExist(end+1) =  0;
    PreFixList{end+1}   = 'PV_pGM';                             TemplateNameList{end+1}    = 'PV_pGM';      SessionsExist(end+1) =  0;
    PreFixList{end+1}   = 'PV_pWM';                             TemplateNameList{end+1}    = 'PV_pWM';      SessionsExist(end+1) =  0;
	PreFixList{end+1}   = 'PV_pCSF';                            TemplateNameList{end+1}    = 'PV_pCSF';      SessionsExist(end+1) =  0;
    PreFixList{end+1}   = 'R1';                                 TemplateNameList{end+1}    = 'R1';          SessionsExist(end+1) =  0;
    %
    % ASL images
    % SessionsExist should be 0 for structural scans, since we only use 1
    % per patient, irrespective of number of functional scans

    PreFixList{end+1}   = ['q' x.P.CBF];                        TemplateNameList{end+1}  = 'CBF';                           SessionsExist(end+1)    = 1;
    PreFixList{end+1}   = 'MaskVascular';                       TemplateNameList{end+1}  = 'MaskVascular';                  SessionsExist(end+1)    = 1;
    PreFixList{end+1}   = 'rMaskSusceptibility';                TemplateNameList{end+1}  = 'MaskSusceptibility';            SessionsExist(end+1)    = 1;
    PreFixList{end+1}   = 'SliceGradient';                      TemplateNameList{end+1}  = 'SliceGradient';                 SessionsExist(end+1)    = 1;
    PreFixList{end+1}   = ['q' x.P.CBF '_masked'];              TemplateNameList{end+1}  = 'CBF_masked';                    SessionsExist(end+1)    = 1;

    PreFixList{end+1}   = 'mean_control_beforeMoCo';            TemplateNameList{end+1}  = 'mean_control_beforeMoCo';       SessionsExist(end+1)    = 1;
    PreFixList{end+1}   = 'SD_control_beforeMoCo';              TemplateNameList{end+1}  = 'SD_control_beforeMoCo';         SessionsExist(end+1)    = 1;
    PreFixList{end+1}   = 'SNR_control_beforeMoCo';             TemplateNameList{end+1}  = 'SNR_control_beforeMoCo';        SessionsExist(end+1)    = 1;
    PreFixList{end+1}   = 'mean_PWI_beforeMoCo';                TemplateNameList{end+1}  = 'mean_PWI_beforeMoCo';           SessionsExist(end+1)    = 1;
    PreFixList{end+1}   = 'SD_PWI_beforeMoCo';                  TemplateNameList{end+1}  = 'SD_PWI_beforeMoCo';             SessionsExist(end+1)    = 1;
    PreFixList{end+1}   = 'SNR_PWI_beforeMoCo';                 TemplateNameList{end+1}  = 'SNR_PWI_beforeMoCo';            SessionsExist(end+1)    = 1;

    PreFixList{end+1}   = 'PWI';                                TemplateNameList{end+1}  = 'PWI';              SessionsExist(end+1)    = 1;
    PreFixList{end+1}   = 'M0';                                 TemplateNameList{end+1}  = 'M0';               SessionsExist(end+1)    = 1;
    PreFixList{end+1}   = 'noSmooth_M0';                        TemplateNameList{end+1}  = 'noSmooth_M0';      SessionsExist(end+1)    = 1;
    PreFixList{end+1}   = 'mean_control';                       TemplateNameList{end+1}  = 'mean_control';     SessionsExist(end+1)    = 1;
    PreFixList{end+1}   = 'SD_control';                         TemplateNameList{end+1}  = 'SD_control';       SessionsExist(end+1)    = 1;
    PreFixList{end+1}   = 'SNR_control';                        TemplateNameList{end+1}  = 'SNR_control';      SessionsExist(end+1)    = 1;

    PreFixList{end+1}   = 'SD';                                 TemplateNameList{end+1}  = 'SD';               SessionsExist(end+1)    = 1;
    PreFixList{end+1}   = 'SNR';                                TemplateNameList{end+1}  = 'SNR';              SessionsExist(end+1)    = 1;
    PreFixList{end+1}   = 'TT';                                 TemplateNameList{end+1}  = 'TT';               SessionsExist(end+1)    = 0;
	PreFixList{end+1}   = 'ATT';                                TemplateNameList{end+1}  = 'ATT';              SessionsExist(end+1)    = 1;
    PreFixList{end+1}   = 'FoV';                                TemplateNameList{end+1}  = 'FoV';              SessionsExist(end+1)    = 1;
    
    PreFixList{end+1}   = '4V_MAP';                             TemplateNameList{end+1}  = '4V_MAP';           SessionsExist(end+1)    = 0;
    PreFixList{end+1}   = '4V';                                 TemplateNameList{end+1}  = '4V';               SessionsExist(end+1)    = 0;
    PreFixList{end+1}   = 'CoW_MAP';                            TemplateNameList{end+1}  = 'CoW_MAP';          SessionsExist(end+1)    = 0;
    PreFixList{end+1}   = 'HEMI';                               TemplateNameList{end+1}  = 'HEMI';             SessionsExist(end+1)    = 0;    
    PreFixList{end+1}   = 'TASL';                               TemplateNameList{end+1}  = 'TASL';             SessionsExist(end+1)    = 0;        
	PreFixList{end+1}   = 'Tex';                                TemplateNameList{end+1}  = 'Tex';              SessionsExist(end+1)    = 1;
	PreFixList{end+1}   = 'ABV';                                TemplateNameList{end+1}  = 'ABV';              SessionsExist(end+1)    = 1;
    % % PM: Let this search for different scantypes in /PopDir NIfTIs, & run within those
end


% ----------------------------------------------------------------------------------------------------
%% 2. Iterate over scan types & sessions
for iScanType=1:length(PreFixList)
    UnAvailable = 0;

	% Define sessions if relevant for the datatype
	if SessionsExist(iScanType)
		% We are looking for all ASL_X sessions
		x.S.InputDataStr = PreFixList{iScanType};
		[nSessions, ~, listSessions] = xASL_adm_GetPopulationSessions(x, false); % run this without verbosity
	else
		% For this scantype, there are not ASL_X sessions, so we set a dummy single session
		listSessions = {'ASL_1'};
		nSessions = 1;
	end
    
    fprintf('%s\n', ['Searching ' TemplateNameList{iScanType} ' images:']);
    for iSession=1:nSessions % iterate over sessions

        if iSession==1 && ~SessionsExist(iScanType)
                % For structural scans, there is no session appendix
                SessionAppendix = '';
                bProceedThisSession = 1;
        elseif iSession>1  && ~SessionsExist(iScanType)
                % For structural scans, there are no sessions>1, so
                % skip this
                bProceedThisSession = 0;
        elseif SessionsExist(iScanType)
                SessionAppendix = ['_' listSessions{iSession}];
                bProceedThisSession = 1;
        end

        if bProceedThisSession

            % ----------------------------------------------------------------------------------------------------
            % Predefine & clear memory
            IM = {0};
            IM2noMask = {0};
            LoadFiles{1} = '';
            LoadFiles{2} = '';
            UnAvailable = 0;
            NoImageN = 1;
            % Searching for available images
            
            if size(x.S.SetsID, 1) ~= x.dataset.nSubjects * x.dataset.nSessions
                error('Mismatch between x.S.SetsID & x.dataset.nSubjects * x.dataset.nSessions');
            end
            
            LoadSetsID = false(size(x.S.SetsID, 1), 1);
            
            AnyBilateralFound = false;
            AnyUnilateralFound = false;
            
            % ----------------------------------------------------------------------------------------------------
            %% 3. Check availability images          
            for iSubject = 1:x.dataset.nSubjects
                SubjSess = (iSubject-1)*nSessions + iSession;
                xASL_TrackProgress(SubjSess, x.dataset.nSubjects* nSessions);
                PathNII = fullfile(x.D.PopDir,[PreFixList{iScanType} '_' x.SUBJECTS{iSubject} SessionAppendix '.nii']);
                PathNII_Left = fullfile(x.D.PopDir,[PreFixList{iScanType} '-L_' x.SUBJECTS{iSubject} SessionAppendix '.nii']);
                PathNII_Right = fullfile(x.D.PopDir,[PreFixList{iScanType} '-R_' x.SUBJECTS{iSubject} SessionAppendix '.nii']);

                % Track if bilateral maps exist
                if xASL_exist(PathNII, 'file')
                    AnyBilateralFound = true; % use this across all subjects/sessions
                    ExistBilateral = true; % use this one here
                else
                    ExistBilateral = false;
                end
                % Track if unilateral (both left & right) maps exist
                if xASL_exist(PathNII_Left, 'file') && xASL_exist(PathNII_Right, 'file')
                    AnyUnilateralFound = true; % use this across all subjects/sessions
                    ExistUnilateral = true; % use this one here
                else
                    ExistUnilateral = false;
                end
                
                if ExistBilateral
                    % If exist, add this subject/image to the list
                    LoadFiles{1}{end+1, 1} = PathNII;
                    LoadSetsID(SubjSess, 1) = true;                         
                elseif ExistUnilateral
                    % same here
                    LoadFiles{1}{end+1, 1} = PathNII_Left;
                    LoadFiles{2}{end+1, 1} = PathNII_Right;
                    LoadSetsID(SubjSess, 1) = true;               
                else
                    % if doesnt exist, dont add to the list
                    UnAvailable = UnAvailable+1;
                    NoImageN = NoImageN+1;
                end
            end
            
            if isempty(LoadFiles{1})
                fprintf('\n%s',['No ' PreFixList{iScanType} ' NIfTIs found, skipping...']);
            elseif AnyBilateralFound && AnyUnilateralFound
                warning('\n%s',['Both bilateral & unilateral ' PreFixList{iScanType} ' NIfTIs found, please remove one of these, skipping...']);
            else
                %% 4. Load images
                if AnyUnilateralFound
                    fprintf('%s\n', 'Unilateral (left and right) images detected');
                    LoadString = {'left' 'right'};
                else
                    LoadString = {'bilateral'};
                end

                % determine whether we load one image per subject or one
                % image per session (== multiple per subject)
                if SessionsExist(iScanType)
                    nSize = x.dataset.nSubjects * x.dataset.nSessions;
                else
                    nSize = x.dataset.nSubjects;
                end

                if bSkipWhenMissingScans && UnAvailable>0.10*nSize % we can allow for 10% unavailable scans
                    fprintf('\n%s',['More than 10% missing ' PreFixList{iScanType} ' files, skipping...']);
                else

                    % If we load both left & right (unilateral) images, IM & IM2noMask
                    % becomes 2 cells, if we load only bilateral images, IM
                    % & IM2noMask have 1 cell
                    bProceedComputationMaps = 1;
                    
                    for iCell=1:2
                        if iCell==2 && AnyBilateralFound
                            % we skip the second cell for bilateral images
                        else
                            fprintf('%s', ['Loading ' LoadString{iCell} ' images:   ']);
                            nLoad = size(LoadFiles{iCell}, 1);
                            IM{iCell} = zeros(Size1, nLoad, 'single'); % pre-allocating for speed
                            if bSaveUnmasked; IM2noMask{iCell} = zeros(121,145,121,nLoad, 'single'); end

                            for iLoad=1:nLoad % add images
                                xASL_TrackProgress(iLoad, nLoad);
                                tempIM = xASL_io_Nifti2Im(LoadFiles{iCell}{iLoad, 1});
                                tempImColumn = xASL_im_IM2Column(tempIM, x.S.masks.WBmask);

                                if iLoad>1 && ~(size(tempImColumn,1)==size(IM{iCell},1))
                                    warning(['Wrong size:' LoadFiles{iCell}{iLoad,1}]);
                                    bProceedComputationMaps = 0; % proceed with next ScanType
                                else % add the image
                                    IM{iCell}(:,iLoad) = tempImColumn;
                                    if bSaveUnmasked; IM2noMask{iCell}(:,:,:,iLoad) = tempIM; end
                                end
                            end
                            fprintf('\n');
                            
                            % clip below zero for visualization
                            IM{iCell}(IM{iCell}<0) = 0;
                            if bSaveUnmasked; IM2noMask{iCell}(IM2noMask{iCell}<0) = 0; end                            
                        end
                    end
                    
                    CurrentSetsID = x.S.SetsID(LoadSetsID, :);
                    
                    if bProceedComputationMaps
                        % initialize image indices that will be included
                        NotOutliers = true(1, size(IM{1}, 2));
                        
                        % ----------------------------------------------------------------------------------------------------
                        %% 5. Remove outliers
                        if bRemoveOutliers
                            % Exclude outliers
                            for iCell=1:length(IM)
                                NotOutliersThisCell = xASL_stat_RobustMean(IM{iCell})';
                                NotOutliers = NotOutliers & NotOutliersThisCell;
                            end
                        end
                        TempOutliers = 1:size(IM{1}, 2);
                        NotOutliers = TempOutliers(NotOutliers);

                        NameIM = [TemplateNameList{iScanType} x.S.TemplateNumberName];
                        
                        % ----------------------------------------------------------------------------------------------------
                        %% 6. Compute templates for all subjects together (only for bilateral images)
                        if length(IM)==1
                            xASL_wrp_CreatePopulationTemplates_Computation(IM{1}(:, NotOutliers), NameIM, x, FunctionsAre, true, SmoothingFWHM);

                            if bSaveUnmasked
                                xASL_wrp_CreatePopulationTemplates_Computation(IM2noMask{1}(:,:,:, NotOutliers), NameIM, x, FunctionsAre, false, SmoothingFWHM);
                            end
                        end
                        
                        % ----------------------------------------------------------------------------------------------------
                        %% 7. Compute templates for individual sets
                        if ~bCompute4Sets
                            % not requested, skipping
                        elseif bCompute4Sets && isempty(Sets2Check)
                            fprintf('\n');
                            warning('There are no sets that we can create statistical maps for');
                        else
                            for iSet=1:length(Sets2Check)
                                % First validate that this set doesnt have continuous data
                                if x.S.Sets1_2Sample(Sets2Check(iSet))==3
                                    warning(['Cannot create maps for non-ordinal set ' x.S.SetsName{Sets2Check(iSet)} ', skipping'])
                                else
                                    % run an iteration for a subset
                                    xASL_wrp_CreatePopulationTemplates4Sets(x, bSaveUnmasked, bRemoveOutliers, FunctionsAre, Sets2Check(iSet), IM, IM2noMask, iScanType, listSessions, SessionsExist, iSession, TemplateNameList, CurrentSetsID, SmoothingFWHM);
                                end
                            end % iSet=1:length(Sets2Check)
                        end % if bComputeSets
                    end % if bProceedComputationMaps
                end % bSkipWhenMissingScans && UnAvailable>0.10*nSize
            end % bSkipWhenMissingScans && isempty(LoadFiles)
        end % if bProceedThisSession
    end % for iSession=1:nSessions
    fprintf('\n');
    if UnAvailable>0
        fprintf('%s\n',[num2str(UnAvailable) ' ' PreFixList{iScanType} ' files missing']);
    end
end % for iScanType=1:length(PreFixList)


end






%% ===================================================================================
%% ===================================================================================
function xASL_wrp_CreatePopulationTemplates4Sets(x, bSaveUnmasked, bRemoveOutliers, FunctionsAre, Set2Check, IM, IM2noMask, iScanType, listSessions, SessionsExist, iSession, TemplateNameList, CurrentSetsID, SmoothingFWHM)
%xASL_wrp_CreatePopulationTemplates4Sets Subfunction that creates the parametric images for subsets
%
% INPUT:
%   x                 - structure containing fields with all information required to run the population module (REQUIRED)
%   bSaveUnmasked     - allows saving the same images without masking (OPTIONAL, DEFAULT=true)
%   bRemoveOutliers   - This parameter makes robust statistics, by removing
%                       outliers (OPTIONAL, DEFAULT=false)
%   FunctionsAre      - two cells, with functions, and with names
%                       (OPTIONAL, DEFAULT is mean and SD). Same setup as
%                       "SpecificScantype" above). First should be a
%                       functionhandle, e.g. @xASL_stat_MeanNan (not a
%                       string)
%   Set2Check         - the set for which parametric maps will be computed 
%   IM                - Cell structure with either one (bilateral images) or two
%                       (unilateral left or right hemispheres separately).
%                       Each cell contains the image matrix, masked/compressed into columns, 
%                       with [IntraMask N] (N==subjects/sessions)
%   IM2noMask         - same as IM but without the masking/compression, so 4D
%                       matrices with [X Y Z N]
%   iSet              - index of current set
%   iScanType         - index of the scan types for which parametric maps are
%                       computed
%   listSessions      - List of the sessions to use
%   SessionsExist     - boolean indicating if sessions exist for this
%                       ScanType
%   iSession          - index of session for which parametric maps are computed
%                       (default = 1, also if there are no sessions)
%   TemplateNameList  - output name for the parametric map NIfTI file
%   CurrentSetsID     - table of values for sets that are currently iterated over
%   SmoothingFWHM     - Full-Width-Half-Maximum in [X Y Z] voxels for smoothing of the output image
%                       (OPTIONAL, DEFAULT = [0 0 0] (i.e. no smoothing)%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This subfunction iterates over different sets/categories/options, to
% create parametric maps for each of them. E.g. if SetsOptions{1}{1} =
% 'Cases' and SetsOptions{1}{2} = 'Controls', this subfunction will create
% separate parametric maps for cases and controls. This function needs
% categorical/ordinal sets, it will be skipped for continuous sets (e.g.
% 'age'). Subject/sessions with 'n/a' will be skipped.
%
% Optionally, this subfunction performs a left-right flip, to average
% hemispheres. When for a given image type, sufficex with both '-left' and '-right'
% exist, these will be loaded accordingly as unilateral hemispheres, such
% that the average hemisphere is unilateral.
%
% The following options should be provided for each subject/session:
% in case left & right hemispheres are split between two NIfTIs:
%   'n/a' = don't include this subject/session
%   'left' = include the left hemisphere only for this subject/session
%   'right' = include the right hemisphere only for this subject/session
%   'both' = include both hemispheres for this subject/session
% in case left & right hemispheres are in the same NIfTI image (not split):
%   'n/a' = don't include
%   'left' = include a non flipped image
%   'right' = include a flipped image
%   'both' = both include a flipped and non-flipped image
%
% This subfunction performs the following steps:
% 0. Admin (define sets and input validations)
% 1. For indices with 'both', clone them so they have both left and right
% 2. Flip images with index 'right'
% 3. Reset SetOptions to inclusion/NaN instead of left/right/NaN
% 4. iterate over the options/categories of this set, to create parametric maps
% __________________________________
% Copyright 2015-2023 ExploreASL


% Get CurrentSet
CurrentSet = CurrentSetsID(:,Set2Check);

% ----------------------------------------------------------------------------------------------------
%% 0. Admin (define sets and input validations)

% Obtain unique options for this set
TempUnique = unique(CurrentSet);
UniqueSet = TempUnique(~isnan(TempUnique));
% Now if we have a 0, switch it to 1
% This can happen if we have a set with string options, where
% the first string is a 0, and then translated into a number.
% This should be then index==1 for SetsOptions
IndexZero = find(UniqueSet==0);
if ~isempty(IndexZero) && sum(UniqueSet==1)==0
    UniqueSet(IndexZero) = 1;
end

SetOptions = x.S.SetsOptions{Set2Check};

% if the only options are left & right (and n/a for missing), assume that this is
% a request to flip hemispheres for the 'right' ones
bFlipHemisphere = min(cellfun(@(y) ~isempty(regexpi(y, '^(.|)(left|right|l|r|n/a|nan|both|bothleftright)(.|)$')), SetOptions));

if length(IM)==1
    bUnilateralImages = false;
elseif length(IM)==2 && bFlipHemisphere
    fprintf('%s\n', 'Processing unilateral images (either left or right)');
    bUnilateralImages = true;
    if ~isequal(size(IM{1}), size(IM{2}))
        warning('Inconsistent size left-right images, skipping');
        return;
    end
elseif length(IM)==2 && ~bFlipHemisphere
    warning('Detected unilateral images but left-right designation missing, skipping');
    return;
else
    error('Incorrect IM matrix size, skipping');
end

if bFlipHemisphere
    fprintf('%s\n', 'Hemisphere encoding (left-right designations) detected, flipping images with designation right');
    
    % ----------------------------------------------------------------------------------------------------
    %% 1. For indices with 'both', clone them so they have both left and right
    % First we clone the image to the end as new image, with the
    % designation 'right', and reset the current designation to 'left'
    % (instead of 'both') so they are both averaged in
    IndexLeft = find(cellfun(@(y) ~isempty(regexpi(y, '^(.|)(left|l)(.|)$')), SetOptions));
    IndexRight = find(cellfun(@(y) ~isempty(regexpi(y, '^(.|)(right|r)(.|)$')), SetOptions));    
    IndexBoth = find(cellfun(@(y) ~isempty(regexpi(y, '^(.|)(both|bothleftright)(.|)$')), SetOptions));
    
    if isempty(IndexLeft)
        SetOptions{end+1} = 'left';
        IndexLeft = length(SetOptions);
    end
    if isempty(IndexRight)
        SetOptions{end+1} = 'right';
        IndexRight = length(SetOptions);
    end
    if isempty(IndexBoth)
        SetOptions{end+1} = 'both';
        IndexBoth = length(SetOptions);
    end    
    
    ImagesBoth = find(CurrentSet==IndexBoth);
    % loop over these indices
    if ~isempty(ImagesBoth)
        fprintf('%s\n', 'Detected designation both, using both left and right hemispheres from these images');
        for iBoth=1:length(ImagesBoth)
            % add new index containing right
            CurrentSet(end+1) = IndexRight;
                
            % set current index to left (instead of both)
            CurrentSet(ImagesBoth(iBoth)) = IndexLeft;
            % copy current image to new image index
            for iIM=1:length(IM)
                IM{iIM}(:,end+1) = IM{iIM}(:,ImagesBoth(iBoth));
            end
        end
    end
    
    % ----------------------------------------------------------------------------------------------------
    %% 2. Flip images with index 'right'
    Images2Flip = CurrentSet==IndexRight;
    
    fprintf('Flipping images:   ');
    
    for iImage=1:size(IM{1}, 2) % we iterate this, to avoid using large memory
        xASL_TrackProgress(iImage, size(IM{1},2));
        if Images2Flip(iImage)
            if bUnilateralImages
                tIM = IM{2}(:,iImage); % flip right image (to left)
            else
                tIM = IM{1}(:,iImage); % flip bilateral image (right-left direction)
            end
            tIM = flip(xASL_im_Column2IM(tIM, x.S.masks.WBmask), 1);
            IM{1}(:,iImage) = xASL_im_IM2Column(tIM, x.S.masks.WBmask);

            if bSaveUnmasked && bUnilateralImages
                IM2noMask{1}(:,:,:,iImage) = flip(IM2noMask{2}(:,:,:,iImage), 1); % flip right image (to left)
            elseif bSaveUnmasked && ~bUnilateralImages
                IM2noMask{1}(:,:,:,iImage) = flip(IM2noMask{1}(:,:,:,iImage), 1); % flip bilateral image (right-left direction)
            end
        end
    end
    fprintf('\n');
    
    % ----------------------------------------------------------------------------------------------------
    %% 3. Reset SetOptions to inclusion/NaN instead of left/right/NaN
    % Now all images are flipped to IM{1} & IM2noMask{1}, for both
    % bUnilateralImages (IM has 2 cells) & ~bUnilateralImages (IM has 1 cell)

    % now we change the set options & ID to inclusion instead of left/right
    % to fool the subsequent processing that should split groups
    % instead of splitting left/right, the subsequent processing will
    % include all (as they may have been flipped above)
    IndexInclusion = find(cellfun(@(y) ~isempty(regexpi(y, '^(.|)(left|right|l|r)(.|)$')), SetOptions)); % these are the indices for inclusion (either left or right)
    SetID = ~max(CurrentSet==IndexInclusion, [], 2)+1;
    SetOptions = {'' 'n/a'}; % include ones, exclude twos
    UniqueSet = [1;2];
else
    % Don't apply the hemispheric averaging by flipping, simply keep as is
    SetOptions = lower(x.S.SetsOptions{Set2Check});
    SetID = CurrentSet;
end

% ----------------------------------------------------------------------------------------------------
%% 4. iterate over the options/categories of this set, to create parametric maps
for iU=1:length(UniqueSet)
    HasNaN = ~isempty(regexpi(SetOptions{UniqueSet(iU)}, '(n/a|nan)')); % skipping NaNs, consider them as outside of a group
    CannotHaveSessions = ~SessionsExist(iScanType) && ~isempty(regexpi(x.S.SetsName{Set2Check}, 'session')); % This SET is for sessions, but there are no sessions for this scantype, so skip this

    if ~HasNaN && ~CannotHaveSessions
        try
            WithinGroup = SetID==UniqueSet(iU);

            if ~SessionsExist(iScanType) && length(WithinGroup) == x.dataset.nSubjects * x.dataset.nSessions
                % if no sessions exist, but "WithinGroup" definition was
                % based on all subject/sessions, then correct this
                CurrSess = (1:length(listSessions):x.dataset.nSubjects * x.dataset.nSessions)' + (iSession-1);
                WithinGroup = WithinGroup(CurrSess);
            end

            % select those within the set/group only
            CurrentIM = IM{1}(:,WithinGroup);

            if bRemoveOutliers
                % Exclude outliers
                NotOutliers = find(xASL_stat_RobustMean(CurrentIM))';
            else
                NotOutliers = 1:size(CurrentIM,2);
            end

            % compute maps
            NameIM = [TemplateNameList{iScanType} '_' x.S.SetsName{Set2Check} '_' SetOptions{UniqueSet(iU)} '_n' num2str(length(NotOutliers))];
            xASL_wrp_CreatePopulationTemplates_Computation(CurrentIM(:,NotOutliers), NameIM, x, FunctionsAre, true, SmoothingFWHM);
            if bSaveUnmasked
                CurrentIM2noMask = IM2noMask{1}(:,:,:,WithinGroup);
                xASL_wrp_CreatePopulationTemplates_Computation(CurrentIM2noMask(:,:,:,NotOutliers),NameIM, x, FunctionsAre, false, SmoothingFWHM);
            end
        catch ME
            warning('Getting set didnt work');
            fprintf('%s\n',ME.message);
        end % try
    end % ~HasNaN && ~CannotHaveSessions
end % iU=1:length(UniqueSet)


end







%% ===================================================================================
%% ===================================================================================
function xASL_wrp_CreatePopulationTemplates_Computation(IM, NameIM, x, FunctionsAre, bMask, SmoothingFWHM)
%xASL_wrp_CreatePopulationTemplates_Computation Subfunction that computes the parametric maps
%
% INPUT:
%   IM              - Image matrix, either masked/compressed into columns, with [IntraMask N] (N==subjects/sessions)
%                     or not, with [X Y Z N]
%   NameIM          - name of the scan type/image
%   x               - structure containing fields with all information required to run the population module (REQUIRED)
%   FunctionsAre    - two cells, with functions, and with names
%                     (OPTIONAL, DEFAULT is mean and SD). Same setup as
%                     "SpecificScantype" above). First should be a
%                     functionhandle, e.g. @xASL_stat_MeanNan (not a
%                     string)
%   bMask           - boolean specifying if IM are masked (2D) or unmasked (4D)
%                     as indicated above
%   SmoothingFWHM   - Full-Width-Half-Maximum in [X Y Z] voxels for smoothing of the output image
%                     (OPTIONAL, DEFAULT = [0 0 0] (i.e. no smoothing)%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This subfunction computes the parametric maps, and saves
% them as NIfTI, with the following steps:
% 0. Admin
% 1. Remove empty images (those for which sum(isfinite())==0)
% 2. Specify the desired function (e.g. 'sum', 'mean', 'stdev' etc)
% 3. Compute the map
% 4. Smooth the map (if requested)
% 5. Save the map as NIfTI
% __________________________________
% Copyright 2015-2020 ExploreASL


    % ----------------------------------------------------------------------------------------------------
    %% 0. Admin
    if isempty(IM)
        warning(['No valid images ' NameIM ' found, skipping']);
        return;
    end

    if bMask
        iDim = 2; % image data as column
    else
        iDim = 4; % image data as image
    end

    % ----------------------------------------------------------------------------------------------------
    %% 1. Remove empty images (those for which sum(isfinite())==0)
    UseIM = ones(size(IM,iDim),1);
    for iM=1:size(IM,iDim)
        if bMask
            CheckSum = xASL_stat_SumNan(IM(:,iM))==0;
        else
            CheckSum = xASL_stat_SumNan(xASL_stat_SumNan(xASL_stat_SumNan(IM(:,:,:,iM))))==0;
        end
        if CheckSum
            UseIM(iM) = 0;
        end
    end

    if bMask
        IM = IM(:,logical(UseIM));
    else
        IM = IM(:,:,:,logical(UseIM));
    end

    fprintf(['Computing ' NameIM ' parametric map(s)...\n']);

    for iFunction=1:length(FunctionsAre{1})
        
        % ----------------------------------------------------------------------------------------------------
        %% 2. Specify the desired function (e.g. 'sum', 'mean', 'stdev' etc)
        FunctionHandle = FunctionsAre{1}{iFunction};
        
        % ----------------------------------------------------------------------------------------------------
        %% 3. Compute the map
        if bMask
            ImageIs = xASL_im_Column2IM(FunctionHandle(IM, 2), x.S.masks.WBmask).* x.GradualSkull;
            PathSave = fullfile(x.D.TemplatesStudyDir, [NameIM '_bs-' FunctionsAre{2}{iFunction} '.nii']);
        else
            ImageIs = FunctionHandle(IM, 4);
            PathSave = fullfile(x.D.TemplatesStudyDir, [NameIM '_bs-' FunctionsAre{2}{iFunction} '_Unmasked.nii']);
        end
        
        % ----------------------------------------------------------------------------------------------------
        %% 4. Smooth the map (if requested)
        if max(SmoothingFWHM)>0
            ImageIs = xASL_im_ndnanfilter(ImageIs, 'gauss', SmoothingFWHM, 1); % smooths with about 6 mm FWHM
        end
        
        % ----------------------------------------------------------------------------------------------------
        %% 5. Save the map as NIfTI
        xASL_io_SaveNifti(x.D.ResliceRef, PathSave, ImageIs);
        
        
        % ----------------------------------------------------------------------------------------------------
        %% 6. Save the current masks for QC
        if bMask && iFunction==1 && x.S.bSaveMasks4QC
            % only recommended for checking masks, otherwise this will take
            % up a lot of disk space
            Masks4D = xASL_im_Column2IM(uint8(IM), x.S.masks.WBmask);
            PathSave = fullfile(x.D.TemplatesStudyDir, ['Masks4D_' NameIM '.nii']);
            xASL_io_SaveNifti(x.D.ResliceRef, PathSave, Masks4D);
        end
    end

    
end