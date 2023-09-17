function [RMS] = xASL_adm_CompareDataSets(RefAnalysisRoot,SourceAnalysisRoot,x,type,mutexState)
%xASL_adm_CompareDataSets Compares two ExploreASL datasets for reproducibility
%
% FORMAT:       [RMS] = xASL_adm_CompareDataSets(RefAnalysisRoot,SourceAnalysisRoot,x,type,mutexState)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Compare data sets is used to ...
%
% - type 0: Only save
% - type 1: Save and evaluate
% - type 2: Only evaluate
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     ...
% __________________________________
% Copyright 2015-2020 ExploreASL

%% -----------------------------------------------------------------------
%% Admin
if  exist('x','var')
    bxASL   = true;
else
    bxASL   = false;
end

if nargin < 4 || isempty(type)
	type = 1;
end

if nargin < 5 || isempty(mutexState)
	mutexState = x.mutex.State{end};
end

%% If in ExploreASL, skip reproducibility testing (this function) when not specifically asked
if  bxASL
    if  isfield(x.settings,'bReproTesting')
        if ~x.settings.bReproTesting || isdeployed
            return;
        end
    else
        return;
    end
end

%% For ExploreASL, define reference & source dirs
if     ~exist('RefAnalysisRoot','var') || ~exist('SourceAnalysisRoot','var')
        ExistVars   = 0;
elseif  isempty(RefAnalysisRoot) || isempty(SourceAnalysisRoot)
        ExistVars   = 0;
else
        ExistVars   = 1;
end

if ~ExistVars && bxASL
    RefAnalysisRoot    = fullfile(x.opts.MyPath,'TestDataSet_Processed');
    SourceAnalysisRoot      = x.dir.xASLDerivatives;
    bSaveRMS                = true;
    RMSPath                 = fullfile(x.dir.xASLDerivatives,'RMS_Reproducibility.mat');
end


Root{1}                     = RefAnalysisRoot;
Root{2}                     = SourceAnalysisRoot;

%% Obtain names
RMS{1}                      = SourceAnalysisRoot;

Modules                     = {'Structural' 'ASL' 'Population'};

if bxASL
    RMS{2}                  = x.mutex.ID;
    RMS{3}                  = mutexState;

    CurrentMod              = find(cellfun(@(x) ~isempty(x),regexp(x.mutex.ID(length('xASL_module_')+1:end), Modules)));
    CurrentState            = mutexState(1:3);
else
    CurrentMod = 1;
    CurrentState = '000';
end

%% -----------------------------------------------------------------------
%% Define paths
for iR=1:length(Root) % iR==1 == reference dataset, iR==2 = Source dataset
    Dlist = xASL_adm_GetFileList(Root{iR},'.*', 'List', [0 Inf], true); % find dirs

    % find folders for 1:3 Structural/ASL/population
    SubjI           = find(cellfun(@(x) ~isempty(regexp(lower(x),'sub')), Dlist)); % subject
    PopI            = find(cellfun(@(x) ~isempty(regexp(lower(x),'population')), Dlist)); % population

    ModuleDir{1}{iR}  = fullfile(Root{iR},Dlist{SubjI});  % subject
    ModuleDir{2}{iR}  = fullfile(ModuleDir{1}{iR},'ASL_1'); % ASL
    ModuleDir{3}{iR}  = fullfile(Root{iR},Dlist{PopI});   % population

    for iMod=1:length(ModuleDir) % find NIfTI files
        % This finds NIfTI files that don't start with a number '\D' (those are reserved for the multi-versioning below)
        ModuleList{iMod}{iR} = xASL_adm_GetFileList(ModuleDir{iMod}{iR},'^\D.*\.nii$','FPList',[0 Inf]);
    end
end




%% -----------------------------------------------------------------------
%% Handle multi-version files
% This list of MultiVersion Files will be version-checked: meaning that after each sub-module run,
% a different version is stored. As this is not required for each file, we save space & speed
% by only doing this for the following files.

MultiVersionList            = {'T1' 'FLAIR' 'M0' 'PWI' 'mean_control' 'ASL4D' 'CBF'};


iMod=CurrentMod; % for structural/ASL/population
% replace files by multiversion names
for iR=1:2 % for reference files (iR==1) & source files (iR==2)

    if ~isempty(ModuleList{iMod}{iR})
        % define the multi version files:
        [Fpath, Ffile, Fext]   = cellfun(@(x) xASL_fileparts(x),ModuleList{iMod}{iR},'UniformOutput', false);
        FileIList{iMod}{iR}    = cellfun(@(x) max(strcmp(MultiVersionList,x)),Ffile);

        for iL=1:length(Fpath)
            NewModuleList{iMod}{iR}{iL,1}               = fullfile(Fpath{iL},[CurrentState Ffile{iL} Fext{iL}]);
        end
        NewModuleList{iMod}{iR}(~FileIList{iMod}{iR})   = ModuleList{iMod}{iR}(~FileIList{iMod}{iR}); % keep original names for non-multiversion files
    else
        NewModuleList{iMod}{iR}                         = '';
    end
end



%% -----------------------------------------------------------------------
%% Creating the multi-version files
% fRun this only for the CurrentModule (structural/ASL/population)
% Create multiple versions for or source (iR==2) NIfTIs
% (when IsRef, this is the same as the reference)

OriList2Copy            =    ModuleList{CurrentMod}{2}(FileIList{CurrentMod}{2});
NewList2Copy            = NewModuleList{CurrentMod}{2}(FileIList{CurrentMod}{2});

if type < 2
	% for type 1 and type 0 do saving of new files
	for iL=1:length(NewList2Copy)
		xASL_Copy(OriList2Copy{iL},NewList2Copy{iL},1);

		% Do the same for the .mat sidecar
		[Ppath, Pfile]      = xASL_fileparts(OriList2Copy{iL});
		OriMatFile          = fullfile(Ppath, [Pfile '.mat']);
		[Ppath, Pfile]      = xASL_fileparts(NewList2Copy{iL});
		NewMatFile          = fullfile(Ppath, [Pfile '.mat']);

		if  exist(OriMatFile,'file')
			xASL_Copy(OriMatFile,NewMatFile,1);
		end
	end
end



IsRef  = strcmp(RefAnalysisRoot,SourceAnalysisRoot); % This will be the case when creating the reference dataset

%%  Skip RMS computation & saving if IsRef
if  IsRef
    return;
end



%% -----------------------------------------------------------------------
%% Obtain RMS
if type > 0
	for iMod=1:length(ModuleDir)
		fprintf(['Obtaining RMS for ' Modules{iMod} ' NIfTIs:   ']);
		if  iMod==CurrentMod % compare with the multi-version NIfTIs
			RMS{end+1}                  = xASL_admin_ComputeRMS(NewModuleList{iMod}{1}, NewModuleList{iMod}{2});
		else % compare with the non-version file
			RMS{end+1}                  = xASL_admin_ComputeRMS(ModuleList{iMod}{1}, ModuleList{iMod}{2});
		end
	end

	%% Save RMS
	if bSaveRMS
		if exist(RMSPath,'file')
			RMSload = load(RMSPath, 'RMS');
			RMSload.RMS(end+1,:) = RMS;
			RMS = RMSload.RMS;
		end

		save(RMSPath,'RMS');
	end
end

end







%% -----------------------------------------------------------------------
function [RMS] = xASL_admin_ComputeRMS(RefList,SourceList)
%xASL_admin_ComputeRMS Checks which files have identical names, and obtains their RMS, only for NIfTI files,
% irrespective of their matrix dimensions

if ~isempty(SourceList) && ~isempty(RefList)

    [~, RefFiles, ~]        = cellfun(@(x) xASL_fileparts(x),RefList,'UniformOutput',false);
    [~, SourceFiles,   ~]   = cellfun(@(x) xASL_fileparts(x),SourceList  ,'UniformOutput',false);

    for iN=1:length(SourceList)
        xASL_TrackProgress(iN,length(SourceList));


        %% Print FileName
        RMS{iN,1}           = SourceFiles{iN};

        % First check if the SourceFile exists as RefFile, and which one it is
        RefIndex       = find(cellfun(@(x) strcmp(x,SourceFiles{iN}),RefFiles));

        if      isempty(RefIndex)
                % if file doesn't exist as reference file
                RMS{iN,2}       = 'Reference not existing';
                RMS{iN,3}       = 'Reference not existing';
        elseif ~xASL_exist(RefList{RefIndex}, 'file')
                % if file doesn't exist as template file
                RMS{iN,2}       = 'Reference was expected but missing. Recreate template files (for multiversioning files)';
                RMS{iN,3}       = 'Reference was expected but missing. Recreate template files (for multiversioning files)';
        else % if it does exist
            %% Print image matrix RMS
            RefIm               = xASL_io_Nifti2Im(RefList{RefIndex});
            SourceIm            = xASL_io_Nifti2Im(SourceList{iN});

			[filepath,filename,~] = xASL_fileparts(RefList{RefIndex});

            for iSz=1:9
                SizeRefIm(iSz)          = size(RefIm,iSz);
                SizeSourceIm(iSz)       = size(SourceIm,iSz);
            end

            if ~isequal(SizeRefIm,SizeSourceIm) % if unequal sizes
                RMS{iN,2}   = 'Unequal matrix size';
			else % calculate mean relative voxelwise difference (i.e. AI)

				if strcmp(filename(1:2),'y_')
					% If this is a flow field and we have a cumulative mask, then use that
					MaskFinal = zeros(size(RefIm));
					MaskFinal(ceil(size(MaskFinal,1)/3):ceil(size(MaskFinal,1)/3*2),...
						ceil(size(MaskFinal,2)/3):ceil(size(MaskFinal,2)/3*2),...
						ceil(size(MaskFinal,3)/3):ceil(size(MaskFinal,3)/3*2),:,:) = 1;
					MaskFinal = MaskFinal > 0;
				else
					% Sort template intensities
					SortRefInt              = sort(RefIm(isfinite(RefIm)&(RefIm>0)));
					SortSourceInt           = sort(SourceIm(isfinite(SourceIm)&(SourceIm>0)));
					ThreshRef               = SortRefInt(round(0.3*length(SortRefInt)));
					ThreshSource            = SortSourceInt(round(0.3*length(SortSourceInt)));
					MaskRef                 = RefIm>ThreshRef;
					MaskSource              = SourceIm>ThreshSource;
					MaskFinal               = MaskRef & MaskSource;
				end

                ColumnRef               = RefIm(MaskFinal);
                ColumnSource            = SourceIm(MaskFinal);

                AvIM                    = (ColumnRef+ColumnSource)./2;
                DiffIM                  = ColumnRef-ColumnSource;
                AI                      = abs(DiffIM)./AvIM(:);
				if strcmp(filename(1:2),'y_')
					RMS{iN,2}               = xASL_stat_MedianNan(AI(:)); % percentage
				else
					RMS{iN,2}               = xASL_stat_MedianNan(AI(:)).*100; % percentage
				end
            end

            %% Print orientation matrix RMS
            %  Check first if a .mat orientation matrix sidecar exists
            [Fpath, Ffile, ~]       = xASL_fileparts(RefList{RefIndex});
            RefMatPath              = fullfile(Fpath,[Ffile '.mat']);
            [Fpath, Ffile, ~]       = xASL_fileparts(SourceList{iN});
            SourceMatPath           = fullfile(Fpath,[Ffile '.mat']);

            if  exist(RefMatPath,'file') && exist(SourceMatPath,'file')
                RefNii              = load(RefMatPath,'mat');
                SourceNii           = load(SourceMatPath,'mat');
            else % use the orientation matrix from the NIfTI
                RefNii              = xASL_io_ReadNifti(RefList{RefIndex});
                SourceNii           = xASL_io_ReadNifti(SourceList{iN});
            end

            % Find voxels to compare
			if length(size(MaskFinal))>3
                % puts indices into pointsX pointsY pointsZ coordinates
				[pointsX,pointsY,pointsZ,~] = ind2sub(size(MaskFinal),find(MaskFinal));
			else
				[pointsX,pointsY,pointsZ] = ind2sub(size(MaskFinal),find(MaskFinal));
			end
			% Generate random numbers with the same size as the vector
			iRand = rand(size(pointsX));
			% We want to take every 1/20 of the voxels to calculate the displacement
            % to speed up the process
			iRand = iRand>0.95;

			pointsX = pointsX(iRand);
			pointsY = pointsY(iRand);
			pointsZ = pointsZ(iRand);
			points = [pointsX';pointsY';pointsZ';ones(1,length(pointsX))];

		    % Calculate the new and old position of these points
			DiffOri = 0;
			for ii=1:size(RefNii.mat,3)
				DiffOriSin = RefNii.mat(:,:,ii)*points - SourceNii.mat(:,:,ii)*points;
				% Calculate the mm distance of the shift
				DiffOri(ii) = xASL_stat_MeanNan(sqrt(xASL_stat_SumNan(DiffOriSin(1:3,:).^2,1)));
			end
			RMS{iN,3} = xASL_stat_MeanNan(DiffOri);
        end
    end
else
    fprintf('\b\b\bNo NIfTIs found, skipping, ');
    RMS     = 'No NIfTIs found';
end

fprintf('\n');

end
