function [x] = xASL_stat_PrintStats(x)
%xASL_stat_PrintStats Print overview of (ROI-) data from all subjects/sessions in
%TSV file
%
% FORMAT: [x] = xASL_stat_PrintStats(x)
% 
% INPUT:
%   x                   - struct containing statistical pipeline environment parameters (REQUIRED)
%   x.S.SaveFile        - name of the TSV output file
%   x.S.SetsName        - names of the covariates/sets
%   x.S.NamesROI        - names of the ROIs
%   x.S.SUBJECTID       - names of the subjects (can differ from x.SUBJECTS
%                         if not all subjects had data)
%   x.S.SetsID          - values of the sets
%   x.S.SetsOptions     - options for sets (if the sets are ordinal)
%                         
% OUTPUT:
%   x                   - same as input
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function prints an overview of statistics from
%              data that were acquired per ROI, in a TSV file. It starts by
%              printing covariates (called "Sets"). Rows will be
%              subjects/sessions, columns will be the sets and
%              ROI-statistics.
%              Any missing data will be skipped (setting them to NaN should
%              have happened in a previous function).
%
%              This function performs the following steps:
%
%              1. First remove previous TSV-file, if already existed
%                 printing to a TSV file can be tricky if it is opened by
%                 Excel. Make sure to close previous versions first,
%                 otherwise this part will crash.
%              2. Print overview of sets to TSV
%                 as explained above. Uses subfunction
%                 xASL_stat_CreateLegend to put legends. Aim is to create a
%                 single TSV file that has a proper overview of the data,
%                 & is self-explanatory to those reading/using it.
%              3. Define number of ASL sessions, force to 1 in case of TT or volume metrics
%              4. Print the overview
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: x = xASL_stat_PrintStats(x);
% __________________________________
% Copyright 2015-2020 ExploreASL



%% -----------------------------------------------------------------------------------------------
%% 1) First remove previous TSV-file, if already existed
fclose all; % safety procedure
if isfield(x.S,'summary_fid')
    x.S = rmfield(x.S,'summary_fid'); 
end
    
try
    xASL_delete(x.S.SaveFile);
catch ME
    warning(['Couldnt delete ' x.S.SaveFile ', if it was opened, please close this file first']);
    fprintf('%s\n',['Message: ' ME.message]);
    return;
end




%% -----------------------------------------------------------------------------------------------
%% 2) Print overview of sets to TSV
try
    xASL_adm_CreateDir(fileparts(x.S.SaveFile));
    x.S.FID = fopen(x.S.SaveFile,'wt');
catch ME
    warning(['Couldnt delete ' x.S.SaveFile ', if it was opened, please close this file first']);
    fprintf('%s\n',['Message: ' ME.message]);
    return;    
end

if x.S.FID<0
    fprintf(['x.S.SaveFile is ' x.S.SaveFile '\n']);
    warning('Is something wrong with the path we try to save the file, perhaps it is too long?');
end

% Print header
fprintf(x.S.FID,'%s\t', 'SUBJECT');

if isfield(x.S,'SetsName')
    for iSet=1:length(x.S.SetsName)
        fprintf(x.S.FID,'%s\t', x.S.SetsName{iSet});
    end
end

for ii=1:length(x.S.NamesROI)
    fprintf(x.S.FID,'%s\t', x.S.NamesROI{ii});
end
fprintf(x.S.FID,'\n');

x.S.Legend = xASL_stat_CreateLegend(x);

if isfield(x.S,'Legend')
    for ii=1:length(x.S.Legend)
        fprintf(x.S.FID,'%s\t', x.S.Legend{ii});
    end
    fprintf(x.S.FID,'\n');
end



%% -----------------------------------------------------------------------------------------------
%% 3) Define number of sessions, force to 1 in case of TT or volume metrics
SingleSessions = {'volume' 'TT' 'PV_pGM' 'PV_pWM' 'pGM' 'pWM' 'mrc1T1' 'mrc2T1'};

HasSingleSessionOnly = sum(cellfun(@(y) ~isempty(regexp(y,['^' x.S.output_ID])), SingleSessions));

if HasSingleSessionOnly
    nSessions = 1;
else
    FileList = xASL_adm_GetFileList(x.D.PopDir,'.*ASL_\d*\.nii$');
    [~, Ffile] = cellfun(@(x) xASL_fileparts(x),FileList, 'UniformOutput',false);
    SessionsN = unique(cellfun(@(x) x(end-4:end),Ffile, 'UniformOutput',false));
    nSessions = length(SessionsN);
end


%% -----------------------------------------------------------------------------------------------
%% 4) Print overview
for iSubject=1:x.nSubjects
    for iSession=1:nSessions
        iSubjectSession = (iSubject-1)* nSessions +iSession;
        
        % check first if this SubjectSession has data, otherwise skip &
        % issue a warning
        if length(x.S.SUBJECTID)<iSubjectSession || size(x.S.DAT, 1)<iSubjectSession
            warning(['Missing data, skipping printing data for ' x.SUBJECTS{iSubject} '_' x.SESSIONS{iSession}]);
        else
            
            % print subject name
            fprintf(x.S.FID,'%s\t', x.S.SUBJECTID{iSubjectSession, 1});

            % print values for other covariates
            % here we always have nSubjects*nSessions values
            % so use "SubjectSession"
            if isfield(x.S,'SetsID')
                for iPrint=1:size(x.S.SetsID,2)
                    String2Print = x.S.SetsID(iSubjectSession,iPrint);
                    if length(x.S.SetsOptions{iPrint})>1 % we need options
                        if length(x.S.SetsOptions{iPrint}) >= length(unique(x.S.SetsID(:,iPrint)))-2 % allow for zeros & NaNs
                            if isnumeric(String2Print) && (int16(String2Print) == String2Print) && String2Print>0 && x.S.Sets1_2Sample(iPrint)~=3
                                String2Print = x.S.SetsOptions{iPrint}{String2Print};
                            end
                        end
                    end
                    fprintf(x.S.FID,'%s\t', xASL_num2str(String2Print));
                end
            end

            % This part is different for volume or TT, since there will
            % be only 1 value per subject (this will be done by the above in which nSessions is set to 1
            for iPrint=1:size(x.S.DAT,2) % print actual data
                fprintf(x.S.FID,'%s\t', xASL_num2str(x.S.DAT(iSubjectSession, iPrint)));
            end
            fprintf(x.S.FID,'\n');
        end
    end
end

fclose(x.S.FID);
x.S = rmfield(x.S,'SaveFile');


end






function [Legend] = xASL_stat_CreateLegend(x)
%xASL_stat_CreateLegend Create a row with legends to be printed in the TSV file
%containing ASL analysis stats

    if ~isfield(x,'S')
        x.S = struct;
    end

    Legend{1} = 'StudyID';
    
    if isfield(x.S,'SetsName')
        for iSet=1:length(x.S.SetsName)
            switch x.S.SetsName{iSet}
                case {'LongitudinalTimePoint', 'SubjectNList', 'Site'}
                    Legend{iSet+1} = 'integer';
                case 'AcquisitionTime'
                    Legend{iSet+1} = 'hhmmss';
                case {'GM_vol', 'GM_L', 'WM_vol', 'WM_L', 'CSF_vol', 'CSF_L'}
                    Legend{iSet+1} = 'Liter'; % 1 Liter = 10 cm*10 cm * 10 cm
                case {'WMH_vol'}
                    Legend{iSet+1} = 'mL';
                case 'GM_ICVRatio'
                    Legend{iSet+1} = 'ratio GM/ICV';
                case 'GMWM_ICVRatio'    
                    Legend{iSet+1} = 'ratio (GM+WM)/ICV';
                case {'CBF_spatial_CoV', 'PseudoCBF_spatial_CoV'}
                    Legend{iSet+1} = '(SD/mean)';
                case 'CBF_spatial_CoV_norm'
                    Legend{iSet+1} = 'ratio (spatial CoV/pseudoCBF spatial CoV)';
                case 'Age'
                    Legend{iSet+1} = 'years';
                case 'WMH_count'
                    Legend{iSet+1} = 'n lesions (integer)';
                case 'WMH_vol'
                    Legend{iSet+1} = 'mL'; % 1 mL = 10 mm * 10 mm * 10 mm (so 1000 1 mm^3 voxels)
                case 'MeanMotion'
                    Legend{iSet+1} = 'mm';
                otherwise
                    Legend{iSet+1} = '...';
            end
        end
        
        for iSet=1:length(x.S.SetsName)
            if strcmp(Legend{iSet+1},'...') && ~isempty(strfind(x.S.SetsName{iSet},'count'))
                Legend{iSet+1} = 'integer';
            end
        end
    end
    
    if isfield(x.S,'NamesROI') && isfield(x.S,'SaveFile')
        for iM=1:length(x.S.NamesROI)
            if ~isempty(findstr(x.S.SaveFile,'spatialCoV'))
                Legend{end+1} = '% SD/mean';
            elseif ~isempty(findstr(x.S.SaveFile,'CBF'))
                Legend{end+1} = 'mL/100g/min';
            else
                Legend{end+1} = '...';
            end
        end
    else
        warning('Could not find ROI names or filename to save, tsv header may be invalid');
    end

end
