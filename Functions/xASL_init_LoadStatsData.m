function [x] = xASL_init_LoadStatsData(x)
%LoadStatsData This function loads all data
% that can be used for statistical analyses
% It looks for pre-defined data first (e.g. defined
% within DATA_PAR.m)
% & looks for mat-files in root directory to add

% NB: Names should be unique.
% If adding data that hasn't been defined before (e.g. DATA_PAR.m)
% SetsName & SetsOptions will be mat-filename
% Sets1_2Sample will be 3 (continuous data)

% For missing data, use NaNs rather than 9999. 9999 will mock up a
% parametric analysis whereas NaNs will be treated as missing

% NB: Note that NaN is kept as NaN in the x.S.SetsID table for
% continuous data (e.g. age). But for ordinal data, it is named as a
% category. E.g. if sex contains a NaN, it will have categories: "male",
% "female" and "NaN"

    if  x.nSubjects==0
        fprintf('%s\n','No variables loaded because no subjects found');
    else

        
        
        
        %% ------------------------------------------------------------------------------------------------------------
        %% 1) Add pre-defined groups, if not already added

        if  isfield(x,'group')
            for iG=1:length(x.group)
                VarName = x.group{iG}.name;
                VarContent = x.group{iG}.code;
                VarOptions = x.group{iG}.options;
                VarSample = x.group{iG}.sets_1_2_sample;

                x = AddStatisticalVariable(x, VarName, VarContent, VarOptions, VarSample);
            end
        end


        
        
        
        %% ------------------------------------------------------------------------------------------------------------
        %% 2) Add from mat-files, if not already added. Start with site
        SiteMat = fullfile(x.D.ROOT, 'Site.mat');
       
        FileList = xASL_adm_GetFileList( x.D.ROOT, '.*\.mat$','FPList',[0 Inf]);

        if exist(SiteMat,'file')
            FList{1} = SiteMat;
            FList(2:1+length(FileList),1) = FileList;
        else
            FList = FileList;
        end
        
        for iL=1:length(FList)
            [path, fname, ext]    = fileparts(FList{iL});
            if strcmp(fname,'x') || strcmp(fname,'xASL') || strcmp(fname,'dcm2niiCatchedErrors')
                % Skip this variable, is part of ExploreASL functionality
            else
                VarName             = fname;
                VarOptions          = {fname};
                VarSample           = 2; % (datatype, 0=one-sample, 1=paired, 2=two-sample, 3=continuous
                % Here, we assume that the datatype is two-sample. The script below
                % "AddStatisticalVariable" will automatically change this into a
                % continuous variable if length(unique)>5

                MatContent  = load(FList{iL});
                if isfield(MatContent,fname)
                    VarContent  = MatContent.(fname);
                    x = AddStatisticalVariable( x, VarName, VarContent, VarOptions, VarSample, FList{iL});
                else
                    fprintf('%s\n',['mat-file ' fname '.mat did not contain stats variable equal to variable name & was skipped']);
                end
            end
        end

    end
end





%% ------------------------------------------------------------------------------------------------------------
function [x] = AddStatisticalVariable(x, VarName, VarContent, VarOptions, VarSample, PathMat)
%AddStatisticalVariable Adds variable content as SetsID
% NB if there are missing values but you still want to use this SET, put
% missing values in as NaN. Using 9999 will mock up any parametric
% analysis, NaNs will simply be skipped

% Maximum number of ordinal options. This makes sure that if there are
% more options than this number, and variate was loaded as "ordinal", it
% is set to be "continuous"

MaxNumOrdinalOptions    = 15;



IfProcess        = 1;

if  size(VarContent,2)<2
    fprintf('%s\n',['Variable ' VarName ' skipped, did not contain data or only subject names']);
    IfProcess= 0; % skip this variable, doesn't contain information
end    

% Check if already exists, then remove previous set
if isfield(x.S,'SetsName') && IfProcess
    iSet = find(strcmp(x.S.SetsName,VarName));
    if ~isempty(iSet)
        % remove existing set, to add new one
        AllSets = 1:length(x.S.SetsName);
        AllSets = AllSets([1:iSet-1 iSet+1:end]);
        x.S.SetsName = x.S.SetsName(AllSets);
        x.S.SetsID = x.S.SetsID(:,AllSets);
        x.S.SetsOptions = x.S.SetsOptions(AllSets);
        x.S.Sets1_2Sample = x.S.Sets1_2Sample(AllSets);
    end
end

if IfProcess % if variable is valid to process
    % Find index

    if ~isfield(x.S,'SetsName')
        INDEX = 1;
    else
        INDEX = length(x.S.SetsName)+1;
    end

    if ~iscell(VarContent)
        warning('Subjects and or sessions should be defined as strings within cells for stats variables!');
        fprintf('%s\n', ['Skipping mat-file ' PathMat]);
        return;
    end

    %% If SubjectNames are numbers, convert to strings
    for iNum=1:size(VarContent,1)
        if isnumeric(VarContent{iNum,1})
            VarContent{iNum,1}  = num2str(VarContent{iNum,1});
        end
    end

    %% remove whitespaces
    VarContent(:,1) = strtrim(VarContent(:,1));

    % Determine whether group contains subject-wise or 
    % SubjectSession-wise values
    ContainsSessions    = xASL_adm_FindStrIndex( VarContent,x.SESSIONS{1});
    if ~isempty(ContainsSessions) && min(min(ContainsSessions~=0))
        ContainsSessions    = 1;
        SessionColumn       = 2;
        DataColumn          = 3;
    else
        ContainsSessions    = 0;
        DataColumn          = 2;
    end




    %% ------------------------------------------------------------------------------------------------------------
    %% Get unique list of data options % check for no empty lines
    % Either for subjects only, or when including sessions
    nextN=1;
    fprintf('%s\n',['Checking data for ' VarName '...  ']);
    for iS=1:x.nSubjects
        xASL_TrackProgress(iS,x.nSubjects);

        HasEmpty        = 0;
        IsFound         = 0;
        iGsubj          = find(strcmp(VarContent(:,1),x.SUBJECTS{iS}));

        if  isempty(iGsubj)
            % retry
            for iiG=1:size(VarContent,1)
                if  strcmp(VarContent{iiG,1},x.SUBJECTS{iS})
                    iGsubj  = iiG;
                end
            end
        end

        for iSess=1:x.nSessions


            if  DataColumn==3 && length(iGsubj)>1
                % if there are sessions, account for them
                iG          = find(strcmp(VarContent(iGsubj,2),x.SESSIONS{iSess}));
                iG          = iGsubj(iG);
            else
                iG          = iGsubj;
            end



            if  ~isempty(iG)
                IsFound         = 1;
                iG              = iG(1);
                % if subjects or sessions are mentioned multiple times
                % in the DataFile, use only the first instance                    


                if  isempty(VarContent{iG,DataColumn}) % check first if it is not empty
                    fprintf('%s\n',['Subject #' num2str(iG) ', ' VarContent{iG,1} ' is empty for variable ' VarName]);
                    HasEmpty     = 1;
                else
                    TempDataList{nextN}     = VarContent{iG,DataColumn};

                    if  islogical(TempDataList{nextN}) % correct for logical
                        if  TempDataList{nextN}==0
                            TempDataList{nextN} = 0;
                        else
                            TempDataList{nextN} = 1;
                        end
                    end
                    nextN                   = nextN+1;
                end
            end
            if ~IsFound
                HasEmpty    = 1;
            end
        end
    end

    fprintf('\n');

    if  ~exist('TempDataList','var')
        fprintf('%s\n',['Variable ' VarName ' not included because all subjects were missing']);
    else

        % If data is a combination of strings & numbers, convert all to
        % strings
        CountString     = 0;
        CountNumber     = 0;

        for iN=1:length(TempDataList)
                if      isstr(TempDataList{iN})
                        CountString     = CountString+1;
                elseif  isnumeric(TempDataList{iN})
                        CountNumber     = CountNumber+1;
                end
        end

        if  CountString>0 && CountNumber>0 % if both strings & number exists
            % convert all to strings
            for iN=1:length(TempDataList)
                TempDataList{iN}    = num2str(TempDataList{iN});
            end
        end

        try     DataOptionList              = unique( TempDataList );
        catch
                % this workaround is here because "unique" doesnt accept
                % cells with numbers, only numbers or cells with strings
                for iD=1:length(TempDataList)
                    NumericDataList(iD)     = TempDataList{iD};
                end
                DataOptionList              = unique(NumericDataList);
        end

        DataOptionList                      = sort(DataOptionList);




        %% ------------------------------------------------------------------------------------------------------------
        %% Correct DataOptionList for multiple NaNs, if data is numeric
        if  isnumeric(DataOptionList)
            NextN       = 1;
            NaNindex    = 1;
            for iD=1:length(DataOptionList)
                if  isnan(DataOptionList(iD))
                    NaNindex(NextN)     = iD;
                    NextN               = NextN+1;
                end
            end

            if  sum(NaNindex>1)
                DataOptionList          = DataOptionList(1:min(NaNindex));
            end        
        end





        %% ------------------------------------------------------------------------------------------------------------
        %% If data are numbers hidden in strings, convert to numbers
        %% Else, convert strings to ordinal numbers
        for iNum=1:size(VarContent,1)
            if  islogical(VarContent{iNum,DataColumn})
                if  VarContent{iNum,DataColumn}==0
                    VarContent{iNum,DataColumn} = 0;
                else
                    VarContent{iNum,DataColumn} = 1;
                end
            elseif ~isnumeric(VarContent{iNum,DataColumn})
                try
                    VarDataIs       = str2num(VarContent{iNum,DataColumn});
                catch
                    error(['Covariant data of ' VarName ' were not numbers!']);
                end

                if  isempty(VarDataIs)
                    % Assume that data are strings, so take unique strings
                    % and assign numbers to them
                    VarDataIs   = find(strcmp(DataOptionList,VarContent{iNum,DataColumn}));
                 end

                try
                    VarContent{iNum,DataColumn}  = VarDataIs;
                catch
                    error(['Covariant data of ' VarName ' were not numbers!']);
                end                
            end
        end

        DataIsContinuous = 0;
        ListContinuousData = {'AcquisitionTime' 'age' 'MeanMotion' 'GM_vol' 'WM_vol' 'CSF_vol' 'GMWM_ICVRatio' 'GM_ICVRatio' 'WMH_count' 'WMH_vol'...
            'CBF_spatial_CoV' 'CBF_spatial_CoV_norm' 'PseudoCBF_spatial_CoV' 'Yrs_AAO' 'Hematocrit' 'SliceReadoutTime' 'AcquisitionTime' 'Inititial_PLD'...
            'LabelingEfficiency' 'qnt_ATT' 'qnt_T1a' 'qnt_lab_eff'...
            'Blood_T1art' 'GM','deep_WM','WholeBrain','L-ICA','R-ICA','POS','Caudate','Cerebellum','Frontal','Insula','Occipital','Parietal','Putamen'...
            'Temporal','Thalamus','ACA_1','ACA_2','ACA_3','MCA_1','MCA_2','MCA_3','PCA_1','PCA_2','PCA_3','ACA_bilat','MCA_bilat','PCA_bilat','WMH','NAWM'};

         % this is a list of known variables that will never be ordinal, always continuous
        % Check if current VarName is one of these
        DataIsQ = find(strcmp(lower(ListContinuousData),lower(VarName)));
        if ~isempty(DataIsQ)
            DataIsContinuous = 1;
            VarSample = 3; % continuous
        end




        %% ------------------------------------------------------------------------------------------------------------
        %% If data are considered ordinal

        if  length(DataOptionList)<MaxNumOrdinalOptions && ~isempty(DataOptionList) && ~DataIsContinuous
            %% Specify variable options, only if the variate is ordinal

            if  isnumeric(DataOptionList)
                for iD=1:length(DataOptionList)
                    VarOptions{iD}      = num2str(DataOptionList(iD));
                end
            else
                VarOptions              = DataOptionList;
            end

            %% If data are categorical numbers, but still not ordinal numbers, convert to ordinal numbers
            %% This is required because only integer numbers>0 can be indices
            %% This will also avoid zeros (by adding 1), since 0's are no valuable index
            if  isnumeric(DataOptionList)
                if  ~min(DataOptionList==[1:1:length(DataOptionList)]) % if numbers are not continuous integers (== ordinal sequence)
                    for iNum=1:size(VarContent,1)
                        for iD=1:length(DataOptionList)
                            if  VarContent{iNum,DataColumn}==DataOptionList(iD)
                                % this doesn't work for NaNs, hence the
                                % elseif statement
                                VarDataIs = iD;
                            elseif isempty(VarContent{iNum,DataColumn})
                                % do nothing, this will stay empty
                            elseif isnan(VarContent{iNum,DataColumn}) && isnan(DataOptionList(iD))
                                % the NaN will become an option that gets
                                % a unique ID, just as any other DataOption
                                VarDataIs = iD;
                            end
                        end
                        if ~exist('VarDataIs', 'var')
                            VarContent{iNum,DataColumn}         = [];
                        else
                            VarContent{iNum,DataColumn}         = VarDataIs;
                        end
                    end
                end
            end
        end




        %% ------------------------------------------------------------------------------------------------------------
        %% Check whether data is complete for all subjects
        CountAbsent         = 0;

        for iSubj=1:x.nSubjects
            for iSess=1:x.nSessions

                % ID (which name, group etc), all for identification
                iSubjSess                 = (iSubj-1)* x.nSessions +iSess;

                % For each group, create an additional set
                % NB: at least all subjects/sessions that are used should be 
                % present, otherwise this will give an error!
                % User is responsible for imputing/interpolating missing values

                % 1st column should be subject-id column
                % 2nd column should be parameter values
                % or 2nd column is session-id (e.g. 'ASL_1', then
                % 3rd column contains parameter values


                NewVariable(iSubjSess,1)            = NaN; % pre-filling for missing values
                % SUBJECT INDEX
                index1                              = xASL_adm_FindStrIndex( VarContent(:,1), x.SUBJECTS{iSubj});

                if  max(isempty(index1) | index1==0)     % if this subject couldn't be found in variable data
                    AbsentSubjects{CountAbsent+1,1}     = x.SUBJECTS{iSubj};
                    CountAbsent                         = CountAbsent+1;                    

                else
                    if ~ContainsSessions            % If there is no session column
                        for iN=1:length(index1) % this debugs any erroneous double entries
                            VC                          = VarContent{ index1(iN),DataColumn};
                        end
                        NewVariable(iSubjSess,1)        = VC(1);
                    else % if there is a session column

                        VarContent(:,2)   = strtrim(VarContent(:,2)); % remove whitespaces

                        % This code assumes 1 column for x.P.SubjectID & 1 column for x.P.SessionID SESSION INDEX
                        index2                                  = xASL_adm_FindStrIndex( VarContent(index1(:,1),SessionColumn), x.SESSIONS{iSess} );

                        if  isempty(index2) | index2==0 % no session data found fpr this subject/session
                            AbsentSubjects{CountAbsent+1,1}     = x.SUBJECTS{iSubj};
                            CountAbsent                         = CountAbsent+1;
                        elseif  isnumeric(VarContent{ index1(index2),DataColumn}) && ~isempty(VarContent{ index1(index2),DataColumn})  % include this session
                                nIndex                          = index1(index2);
                                for iN=1:length(nIndex) % this debugs any erroneous double entries
                                    VC                          = VarContent{ nIndex(iN),DataColumn};
                                end
                                NewVariable(iSubjSess,1)        = VC(1);
                        else
                            AbsentSubjects{CountAbsent+1,1}     = x.SUBJECTS{iSubj};
                            CountAbsent                         = CountAbsent+1;
                        end % if isempty(index2) | index2==0
                    end % if ~ContainsSessions
                end % if  max(isempty(index1) | index1==0
            end % for iSess=1:x.nSessions
        end % for iSubj=1:x.nSubjects





        %% ------------------------------------------------------------------------------------------------------------
        %% If data is complete, include it in x.S.SETS
        %  Incomplete data is fine, when it is continuous data, to be
        %  filled by NaNs
        NewVariableUnique   = unique(NewVariable(~isnan(NewVariable))); % get unique values, excluding NaNs
        ManyAbsent = CountAbsent>0.05*x.nSubjects || HasEmpty;
        CheckLength = length(NewVariableUnique)<0.25*x.nSubjects;
        CheckLength2 = size(NewVariable, 1)~=size(x.S.SetsID, 1);

        if  (ManyAbsent && CheckLength) || CheckLength2 % don't include this variable because of missing data 
            % Allow for 5% missing data, if data is continuous
            % (assuming continuous data has at least 25% of sample size
            % as unique values
            if exist('AbsentSubjects', 'var')
                fprintf('%s\n',['Variable ' VarName ' not included because the following subjects were missing in this variable:']);
                for iAb=1:size(AbsentSubjects,1)
                    fprintf('%s',[AbsentSubjects{iAb,1} ', ']);
                    if  (iAb/12)==ceil(iAb/12) % 12 subjects per line
                        fprintf('\n'); % next line
                    end
                end
            else
                fprintf('%s\n',['Variable ' VarName ' not included because some data was missing, strangely. Please check']);
            end
            fprintf('\n');
        elseif  size(NewVariable, 1)==size(x.S.SetsID, 1)
            % include variable
            x.S.SetsID(:,INDEX)     = NewVariable;

            % Add variable content, name, option & sets_1_2_sample
            x.S.SetsName{INDEX}         = VarName;
            x.S.Sets1_2Sample(INDEX)  = VarSample;         
            x.S.SetsOptions{INDEX}      = VarOptions;

            % if data seems continuous, force Sets1_2Sample==3

            ListContinuousData = {'AcquisitionTime' 'spatial_CoV' 'ICVRatio' 'age' 'MeanMotion' 'GM_vol' 'WM_vol' 'CSF_vol' 'GMWM_ICVRatio' 'GM_ICVRatio' 'WMH_count' 'WMH_vol'...
                'CBF_spatial_CoV' 'CBF_spatial_CoV_norm' 'PseudoCBF_spatial_CoV' 'Yrs_AAO' 'Hematocrit' 'SliceReadoutTime' 'AcquisitionTime' 'Inititial_PLD'...
                'Blood_T1art' 'GM','deep_WM','WholeBrain','L-ICA','R-ICA','POS','Caudate','Cerebellum','Frontal','Insula','Occipital','Parietal','Putamen'...
                'Temporal','Thalamus','ACA_1','ACA_2','ACA_3','MCA_1','MCA_2','MCA_3','PCA_1','PCA_2','PCA_3','ACA_bilat','MCA_bilat','PCA_bilat','WMH','NAWM'};

            % this is a list of known variables that will never be ordinal, always continuous
            % These are common continuous variables generated by ExploreASL, and should be treated as continuous
            % Check if current VarName is one of these
            DataIsContinuous = ~isempty(find(strcmp(lower(ListContinuousData),lower(VarName))));                

            if DataIsContinuous || (length(unique(x.S.SetsID(:,INDEX)))>MaxNumOrdinalOptions && x.S.Sets1_2Sample(INDEX)~=3)
                x.S.Sets1_2Sample(INDEX) = 3;
                fprintf('%s\n',['Set ' x.S.SetsName{INDEX} ' is assumed to contain continuous data, not ordinal']);
            end




            %% ------------------------------------------------------------------------------------------------------------
            %% Add ordinal variable options, if they are not yet there
            if  x.S.Sets1_2Sample(INDEX)~=3 && length(x.S.SetsOptions{INDEX})==1
                UniqueData      = unique( x.S.SetsID(:,INDEX) );
                for iOp=1:length( UniqueData )
                    x.S.SetsOptions{INDEX}{iOp}      = num2str(UniqueData(iOp));
                end
            end
        end % if  (CountAbsent>0.05*x.nSubjects || HasEmpty) && NewVariableUnique<0.25*x.nSubjects % don't include this variable because of missing data 
    end % if  ~exist('TempDataList','var')
end % if IfProcess % if variable doesn't already exist
   

end