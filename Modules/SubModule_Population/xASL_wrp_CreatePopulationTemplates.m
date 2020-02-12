function xASL_wrp_CreatePopulationTemplates(x, SaveUnmasked, bCompute4Sets, SpecificScantype, bSkipMissingScans, bRemoveOutliers, FunctionsAre)
%xASL_wrp_CreatePopulationTemplates ExploreASL Population module wrapper,
%creates population parametric images for each ScanType
%
% FORMAT: xASL_wrp_CreatePopulationTemplates(x[, SaveUnmasked, Compute4Sets])
%
% INPUT:
%   x            - structure containing fields with all information required to run the population module (REQUIRED)
%   SaveUnmasked - allows saving the same images without masking (OPTIONAL, DEFAULT=true)
%   bCompute4Sets - creates the same parametric images for subsets of the
%                  data (e.g. cohorts, sites, etc) (OPTIONAL, DEFAULT=false)
%   SpecificScantype - allows providing ScanTypes to create templates for:
%                      should be 3 cells:
%                      SpecificScantype{1} = PreFixList (prefix string of the NIfTI, associated with the ScanType, e.g. 'qCBF', 'PV_PGM', etc)
%                      SpecificScantype{2} = TemplateNameList (output name to give to the NIfTI file of the template)
%                      SpecificScantype{3} = SessionsExist (whether sessions exist, if true it uses ASL_1 ASL_2 ASL_n)
%                      
%                      (OPTIONAL, DEFAULT = use predefined ScanTypes)
%   bSkipMissingScans - This parameter allows to choose if we want to
%                       create templates from incomplete datasets or skip
%                       those templates with missing subjects
%                       (OPTIONAL, DEFAULT=true)
%   bRemoveOutliers   - This parameter makes robust statistics, by removing
%                       outliers (OPTIONAL, DEFAULT=false)
%   FunctionsAre      - two cells, with functions, and with names
%                       (OPTIONAL, DEFAULT is mean and SD). Same setup as
%                       "SpecificScantype" above)
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
% EXAMPLE: xASL_wrp_CreatePopulationTemplates(x);
% EXAMPLE for specific scantypes:
%          xASL_wrp_CreatePopulationTemplates(x, [], [], {{'qCBF'} {'CBF'} 1});
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2019 ExploreASL


%% ----------------------------------------------------------------------------------------------------
%  Admin

if nargin<2 || isempty(SaveUnmasked)
    SaveUnmasked = true;
end
if nargin<3 || isempty(bCompute4Sets)
    bCompute4Sets = false;
end
if nargin<5 || isempty(bSkipMissingScans)
    bSkipMissingScans = true;
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

if x.nSubjectsSessions<6
    % With too small datasets, created templated won't be reliable
    fprintf('%s\n',['Too few images (' num2str(x.nSubjectsSessions) ') for creating templates/parametric images, skipping...']);
    return;
end

Size1 = sum(x.WBmask(:));

x.D.TemplatesStudyDir = fullfile(x.D.PopDir, 'Templates');
xASL_adm_CreateDir(x.D.TemplatesStudyDir);

if bCompute4Sets
    % Define sets to create comparative stats for
    NextN = 1;
    for iSet=1:size(x.S.SetsID,2)
        TempUnique = unique(x.S.SetsID(:,iSet));
        nUniqueValues = length(TempUnique(~isnan(TempUnique)));
        if  nUniqueValues>1 && nUniqueValues<6
            Sets2Check(NextN,1) = iSet;
            NextN = NextN+1;
        end
    end
end


%% ----------------------------------------------------------------------------------------------------
%  Here we specify the images/scantypes
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
    PreFixList{end+1}   = 'PV_pGM';                             TemplateNameList{end+1}    = 'PV_pGM';      SessionsExist(end+1) =  0;
    PreFixList{end+1}   = 'PV_pWM';                             TemplateNameList{end+1}    = 'PV_pWM';      SessionsExist(end+1) =  0;
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
    PreFixList{end+1}   = 'FoV';                                TemplateNameList{end+1}  = 'FoV';              SessionsExist(end+1)    = 1;
    % % PM: Let this search for different scantypes in /PopDir NIfTIs, & run within those
end

%% ----------------------------------------------------------------------------------------------------
%  Algorithms start here

% ----------------------------------------------------------------------------------------------------
% Loading data
for iScanType=1:length(PreFixList)
    UnAvailable = 0;
    NextN = 1;

    fprintf(['Searching ' TemplateNameList{iScanType} ' images:   ']);
    for iSession=1:x.nSessions

        if iSession==1 && ~SessionsExist(iScanType)
                % For structural scans, there is no session appendix
                SessionAppendix = '';
                UseThisSession = 1;
        elseif iSession>1  && ~SessionsExist(iScanType)
                % For structural scans, there are no sessions>1, so
                % skip this
                UseThisSession  = 0;
        elseif SessionsExist(iScanType)
                SessionAppendix = ['_' x.SESSIONS{iSession}];
                UseThisSession = 1;
        end

        % ----------------------------------------------------------------------------------------------------
        % Predefine memory
        clear IM LoadFiles NoIm
        UnAvailable = 0;
        NextN = 1;
        NoImageN = 1;
        % Searching for available images
        if UseThisSession
            for iSubject = 1:x.nSubjects
                SubjSess = (iSubject-1)*x.nSessions + iSession;
                xASL_TrackProgress(SubjSess,x.nSubjects*x.nSessions);
                PathNII = fullfile(x.D.PopDir,[PreFixList{iScanType} '_' x.SUBJECTS{iSubject} SessionAppendix '.nii']);

                if xASL_exist(PathNII,'file')
                    LoadFiles{NextN,1} = PathNII;
                    NextN = NextN+1;
                    xASL_io_ReadNifti(PathNII);
                else
                    UnAvailable = UnAvailable+1;
                    NoIm(NoImageN) = SubjSess;
                    NoImageN = NoImageN+1;
                end
            end
        end

        if bSkipMissingScans
            if ~UseThisSession
                continue;
            elseif UseThisSession && ~exist('LoadFiles','var')
                    fprintf('\n%s',['No ' PreFixList{iScanType} ' files found, skipping...']);
                    continue;
            end
        end

        fprintf(', loading images:   ');
        % load data

        if SessionsExist(iScanType)
            nSize = x.nSubjectsSessions;
        else
            nSize = x.nSubjects;
        end

        nLoad = size(LoadFiles,1);
        if bSkipMissingScans        
            if UnAvailable>0.10*nSize % we can allow for 10% unavailable scans
                fprintf('\n%s',['More than 10% missing ' PreFixList{iScanType} ' files, skipping...']);
                continue;
            end
        end

        IM = zeros(Size1,nSize,'single'); % pre-allocating for speed
        if SaveUnmasked; IM2 = zeros(121,145,121,nSize, 'single'); end

        if exist('NoIm','var')
            IM(:,NoIm) = NaN; % fill the non-existing scans
        end
        for iSubject=1:nLoad % add images
            xASL_TrackProgress(iSubject,nLoad);
            tempIM = xASL_io_Nifti2Im(LoadFiles{iSubject,1});
            tempIM(tempIM<0) = 0; % clipping below zero for visualization
            IM(:,iSubject) = xASL_im_IM2Column(tempIM, x.WBmask);
            if SaveUnmasked; IM2(:,:,:,iSubject) = tempIM; end
        end
        fprintf('\n');

        if bRemoveOutliers
            % Exclude outliers
            NotOutliers = find(xASL_stat_RobustMean(IM))';
        else
            NotOutliers = [1:1:size(IM,2)];
        end

        % create the maps
        NameIM = TemplateNameList{iScanType};
        ComputeParametricImages(IM(:,NotOutliers),NameIM, x, FunctionsAre, true);

        if SaveUnmasked
            ComputeParametricImages(IM2(:,:,:,NotOutliers),NameIM, x, FunctionsAre, false);
        end
        % ----------------------------------------------------------------------------------------------------
        % This part checks for individual sets (e.g. create statistic images for each cohort/session etc)
        if bCompute4Sets
            for iSet=1:length(Sets2Check)
                TempUnique = unique(x.S.SetsID(:,Sets2Check(iSet)));
                UniqueSet = TempUnique(~isnan(TempUnique));
                % Now if we have a 0, switch it to 1
                % This can happen if we have a set with string options, where
                % the first string is a 0, and then translated into a number.
                % This should be then index==1 for SetsOptions
                IndZero = find(UniqueSet==0);
                if ~isempty(IndZero) && sum(UniqueSet==1)==0
                    UniqueSet(IndZero) = 1;
                end
                for iU=1:length(UniqueSet)
                    if ~SessionsExist(iScanType) && ~isempty(findstr(x.S.SetsName{Sets2Check(iSet)},'session'))
                        % This SET is for sessions, but there are no sessions for this scantype
                        % So skip this
                    else
                        try
                            WithinGroup = x.S.SetsID(:,Sets2Check(iSet))==UniqueSet(iU);

                            if ~SessionsExist(iScanType) % if no sessions exist, only take current session here
                                CurrSess = [1:x.nSessions:x.nSubjectsSessions]' + (iSession-1);
                                WithinGroup = WithinGroup(CurrSess);
                            end

                            % select those within the set/group only
                            NameIM = [TemplateNameList{iScanType} '_' x.S.SetsName{Sets2Check(iSet)} '_' x.S.SetsOptions{Sets2Check(iSet)}{UniqueSet(iU)}];
                            tempIM = IM(:,WithinGroup);

                            % compute outliers
                            NotOutliers = xASL_stat_RobustMean(tempIM);
                            % compute maps
                            ComputeParametricImages(tempIM(:,NotOutliers), NameIM, x, FunctionsAre, true);
                            if SaveUnmasked
                                tempIM2 = IM2(:,:,:,WithinGroup);
                                ComputeParametricImages(tempIM2(:,:,:,NotOutliers),NameIM, x, false);
                            end
                        end
                    end
                end
            end
        end
    end % for iSession=1:x.nSessions
    fprintf('\n');
    if UnAvailable>0
        fprintf('%s\n',[num2str(UnAvailable) ' ' PreFixList{iScanType} ' files missing']);
    end
end  % for iScanType=1:length(PreFixList)

end


%% ===================================================================================
function ComputeParametricImages(IM, NameIM, x, FunctionsAre, bMask)
%ComputeParametricImages Subfunction that computes the parametric images

    if bMask
        iDim = 2; % image data as column
    else
        iDim = 4; % image data as image
    end

    % Remove empty images
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

    if isempty(IM)
        warning(['No valid images ' NameIM ' found']);
        if bMask
            IM = zeros(sum(x.WBmask(:)),1);
        else
            IM = zeros(size(x.WBmask));
        end
    end
    fprintf(['Computing & saving ' NameIM ' parametric images, ']);

    for iFunction=1:length(FunctionsAre{1})
        FunctionHandle = FunctionsAre{1}{iFunction};
        if bMask
            ImageIs = xASL_im_Column2IM(FunctionHandle(IM, 2),x.WBmask).* x.GradualSkull;
            PathSave = fullfile(x.D.TemplatesStudyDir, [NameIM '_bs-' FunctionsAre{2}{iFunction} '.nii']);
        else
            ImageIs = FunctionHandle(IM, 4);
            PathSave = fullfile(x.D.TemplatesStudyDir, [NameIM '_bs-' FunctionsAre{2}{iFunction} '_Unmasked.nii']);
        end
        xASL_io_SaveNifti(x.D.ResliceRef, PathSave, ImageIs);
    end

end