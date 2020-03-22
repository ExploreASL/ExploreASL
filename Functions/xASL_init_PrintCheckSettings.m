function x = xASL_init_PrintCheckSettings( x )
%xASL_init_PrintCheckSettings Check whether pre-defined settings existed in DATA_PAR.m
% Prints these on the screen

%% Check obligated x
field_symbol    = {'subject_regexp'};

for iField=1:length(field_symbol)
    if ~isfield(x,field_symbol{iField})
        warning(['x.' field_symbol{iField} ' was not defined in DATA_PAR.m!'])
    end
end

if ~isfield(x,'Q')
    x.Q     = struct;
    warning('x.Q didn''nt exist');
% else
%     field_symbol    = {'Initial_PLD','LabelingDuration','LabelingType','BackGrSupprPulses'};
%     for iField=1:length(field_symbol)
%         if ~isfield(x.Q,field_symbol{iField})
%             warning(['x.Q.' field_symbol{iField} ' was not defined in DATA_PAR.m!'])
%         end
%     end
end

if ~isfield(x,'D')
    warning('x.D didn''nt exist');
else
    field_symbol    = {'ROOT'};
    for iField=1:length(field_symbol)
        if ~isfield(x.D,field_symbol{iField})
            warning(['x.D.' field_symbol{iField} ' was not defined in DATA_PAR.m!'])
        end
    end
end

if ~isfield(x,'Quality') || (x.Quality~=0 && x.Quality~=1)
    x.Quality = 1;
    fprintf('%s\n','Default Quality=1 used (optimal quality)');
end
if ~isfield(x,'DELETETEMP') || (x.DELETETEMP~=0 && x.DELETETEMP~=1)
    x.DELETETEMP = 1;
    fprintf('%s\n','Default x.DELETETEMP=1 used (delete files temporarily used for processing)');
end
if isfield(x,'readout_dim')
    if  strcmp(x.readout_dim,'2D') && ~isfield(x.Q,'SliceReadoutTime')
        warning('x.Q.SliceReadoutTime should be defined for a 2D (multi-slice) sequence!');
    end
end

%% Print settings to check
fprintf('%s\n','-------------------------------------------');
fprintf('%s\n\n','ExploreASL will run with following settings:');
fprintf('%s\n\n',['Root folder = ' x.D.ROOT]);

if x.nWorkers>1
    fprintf(['I am worker ' num2str(x.iWorker) '/' num2str(x.nWorkers) '\n']);
    fprintf('Note that the resulting number of scans mentioned below applies only to this worker\n');
end

fprintf('%s\n',[num2str(x.nTotalSubjects) ' scans - ' num2str(x.nExcluded) ' exclusions, resulting in ' num2str(x.nSubjects) ' scans of: ']);

for iT=1:x.nTimePointsTotal
    fprintf('%s\n',['Longitudinal timePoint ' num2str(iT) ' = ' num2str(x.nTimePointTotalSubjects(iT)) ' scans - ' num2str(x.nTimePointExcluded(iT)) ' exclusions = ' num2str(x.nTimePointSubjects(iT)) ' scans']);
end

fprintf('%s\n',['ASL sessions: ' num2str(x.nSessions)]);

fprintf('\n%s\n','Ancillary data, sets:');
if isfield(x.S,'SetsID')
        fprintf('%s\n',[num2str(size(x.S.SetsID,2)) ' sets are defined for ' num2str(size(x.S.SetsID,1)) ' "SubjectsSessions":']);

        for iSet=1:size(x.S.SetsID,2)
            fprintf(['Set ' num2str(iSet) ' = "' x.S.SetsName{iSet} '" options ']);
            for iOption=1:size(x.S.SetsOptions{iSet},2)
                fprintf(['"' x.S.SetsOptions{iSet}{iOption} '"']);
                if iOption~= size(x.S.SetsOptions{iSet},2)
                    fprintf(' & ');
                end
            end
            if      x.S.Sets1_2Sample(iSet)==1
                    fprintf(', codes for paired data');
            elseif  x.S.Sets1_2Sample(iSet)==2
                    fprintf(', codes for two-sample data');
            elseif  x.S.Sets1_2Sample(iSet)==3
                    fprintf([', continuous variate (with ' num2str(length(unique(x.S.SetsID(:,iSet)))) ' unique values)']);
            end                    

            fprintf('\n');
        end
else    
    fprintf('%s\n','No sets are defined');
end

if ~isfield(x,'M0')
%     warning('M0 option missing!');
else
    fprintf('\n%s\n',['M0 option selected is "' num2str(x.M0) '"']);
end
fprintf('%s\n',['x.DELETETEMP = ' num2str(x.DELETETEMP) ' (delete temporary files)']);
fprintf('%s\n',['x.Quality    = ' num2str(x.Quality) ' (0 = fast try-out; 1 = normal high quality)']);

% Warnings
fprintf('\n%s\n\n','---------------------------------------------');

if ~isempty(regexp(x.subject_regexp, '^(\^|)\.\*(\$|)$'))
    warning('Subject regexp not specific! Check that no wrong folders are included as subjects');
end

fprintf('\n');

end

