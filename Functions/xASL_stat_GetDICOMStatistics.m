function xASL_stat_GetDICOMStatistics(x, ScanType, HasSessions)
%xASL_stat_GetDICOMStatistics Collect DICOM metadata in CSV file to check
%
% FORMAT: xASL_stat_GetDICOMStatistics(x, ScanType, HasSessions)
% 
% INPUT:
%   x           - struct containing pipeline environment parameters (REQUIRED)
%   ScanType    - type of data (e.g. T1, FLAIR, ASL) (REQUIRED)
%   HasSessions - true for multiple ASL sessions (OPTIONAL, DEFAULT=false)
%
% OUTPUT: n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This functions prints DICOM metadata (e.g. parameters used
%              for quantification) and collects them in a single csv.
%              Summarizes this for the total population & checks whether outliers exist.
%              Can be useful to detect software upgrades, where only slight
%              parameter changes can hint on quantification changes
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_stat_GetDICOMStatistics(x, 'ASL', true);
% __________________________________
% Copyright 2016-2019 ExploreASL
% PM: Need to change this script into accepting JSON as well as _parms.mat

%% -----------------------------------------------------------------------------------------------
%% Admin
if nargin<3 || isempty(HasSessions)
    HasSessions = false;
end

matFields = {'RepetitionTime' 'EchoTime' 'NumberOfTemporalPositions' 'MRScaleSlope' 'RescaleSlopeOriginal' 'RescaleIntercept' 'Scanner' 'SliceReadoutTime'};
nFields = length(matFields);
SavePath = fullfile(x.D.DICOMparameterDir,[ScanType '_quantification_parameters.csv']);
SummaryFID = fopen(SavePath,'wt');

% Print header
fprintf(SummaryFID,'%s,', 'SUBJECT,SESSION');
for iField=1:nFields
    fprintf(SummaryFID,'%s,', matFields{iField} );
end
fprintf(SummaryFID,'\n' );

% Load & save individual parameter files
fprintf('%s\n','Loading & saving individual parameter files...  ');

for iSubject=1:x.nSubjects
    xASL_TrackProgress(iSubject,x.nSubjects);
    for iSession=1:x.nSessions
        iSubjSess = (iSubject-1)*x.nSessions+iSession;

        fprintf(SummaryFID,'%s,', x.SUBJECTS{iSubject});
        fprintf(SummaryFID,'%s,', x.SESSIONS{iSession});

        if HasSessions
            mat_load = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, x.SESSIONS{iSession}, [ScanType '_parms.mat']);
        else
            mat_load = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, [ScanType '_parms.mat']);
        end

        if exist(mat_load, 'file')
            MatL = load( mat_load );
            parms = MatL.parms;

            % print all fields for subject_session & load them into x.S.par variabele for later calculations
            for iField=1:nFields
                if isfield(parms,matFields{iField}) && ~isempty(parms.(matFields{iField}))
                    
                    fprintf(SummaryFID,'%s,', num2str( parms.(matFields{iField}) ) );
                    if isnumeric( parms.(matFields{iField}) )
                        TempData = parms.(matFields{iField});
                        if length(TempData)>1
                            warning(['parms.' matFields{iField} ' had multiple values']);
                        end
                        x.S.par(iSubjSess,iField) = min(TempData);
                    else
                        x.S.par(iSubjSess,iField) = NaN; % fill empties with NaNs, not zeros!
                    end
                else % if field doesn't exist for a certain SubjectSession, fill empty cell
                     fprintf(SummaryFID,'%s,', 'NaN' );
                     x.S.par(iSubjSess,iField) = NaN; % fill empties with NaNs, not zeros!
                end
            end

        else % fill empties with NaNs, not zeros!
            x.S.par(iSubjSess,1:nFields) = NaN;
        end
        fprintf(SummaryFID,'\n');
    end
end

% Print summary 

fprintf(SummaryFID,'\n\n');
fprintf('\n');

if ~isfield(x.S,'par')
    fprintf('%s\n',['Checking DICOM header values skipped for ' ScanType ' because no images existed']);
else
    fprintf('%s\n',['Printing summary of ' ScanType ' DICOM values & checking for outliers']);

    for iField=1:size(x.S.par,2)
        temp_data = x.S.par(:,iField);
        statField{1}(iField) = xASL_stat_MeanNan( temp_data );
        statField{2}(iField) = xASL_stat_StdNan( temp_data);
        statField{3}(iField) = 100 * statField{2}(iField) / statField{1}(iField) ;
        statField{4}(iField) = statField{1}(iField) + (3 * statField{2}(iField) );
        statField{5}(iField) = statField{1}(iField) - (3 * statField{2}(iField) );
        clear temp_data
    end

    fieldStats = { 'mean' 'SD' '"CV (coeff var, %)"' 'up_threshold (mean+3SD)' 'lo_threshold (mean-3SD)' };

    for iStat=1:length(fieldStats)
        fprintf(SummaryFID,'%s,,', fieldStats{iStat}  );
        for iField=1:size(x.S.par,2)
            fprintf(SummaryFID,'%s,', num2str(statField{iStat}(iField)) );
        end
        fprintf(SummaryFID,'\n' );
    end
end
fclose(SummaryFID);

clear SummaryFID SavePath




%% -----------------------------------------------------------------------------------------------
%% Print outliers

if isfield(x.S,'par')    
    fprintf('%s\n','Printing outliers...  ');
     for iField=1:size(x.S.par,2)
         xASL_TrackProgress(iField,size(x.S.par,2));

        OutlierIndices{iField} = find( x.S.par(:,iField) > statField{4}(iField) | x.S.par(:,iField) < statField{5}(iField) );

        if ~isempty(OutlierIndices{iField})
            % open new csv file
            SavePath = fullfile(x.D.DICOMparameterDir,[ScanType '_outliers_' matFields{iField} '.csv']);
            SummaryFID = fopen(SavePath,'wt');

            % Print header
            fprintf(SummaryFID,'%s\n', ['SUBJECT,SESSION,' matFields{iField}]);

            for ii=1:length( OutlierIndices{iField} )
                iSubject = ceil(OutlierIndices{iField}(ii)/x.nSessions);
                iSession = OutlierIndices{iField}(ii) - ((iSubject-1)*x.nSessions);

                fprintf(SummaryFID,'%s,%s,%s', num2str(x.SUBJECTS{iSubject}),num2str(x.SESSIONS{iSession}),num2str(x.S.par( OutlierIndices{iField}(ii),iField) ) );
                fprintf(SummaryFID,'\n' );
            end

            fclose(SummaryFID);
        end
     end
end

fprintf('\n');

end

