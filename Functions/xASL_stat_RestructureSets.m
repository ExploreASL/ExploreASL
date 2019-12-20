function DataOut = xASL_stat_RestructureSets( DataIn, IDColumn )
%xASL_stat_RestructureSets Restructures a 2D data-matrix into multiple cells with populations
% Data-matrix should be a row for each subject
% Columns represent data points for multiple variables
% Code list & data list should have exact row reference
% Code list should be [first column] numbers of subjects and [second column] code identification.
% Identification code can be any nonzero number, but low numbers are easiest since this function
% restructures populations in cells dimension. e.g. 1 and 2
%
% if a subject's code is missing, put in a NaN, this skips this subject

% PM: in the case that there are groups 1 2 3 and 4, but the data that
% needs to be split doesn't contain group 3, than this will crash.
% So currently all codes need to be present, for the code to run


    codes   = unique(IDColumn(isfinite(IDColumn))); % Get different codes
    nCodes  = length(codes);

    for iCode=1:nCodes
        NextN(iCode)=1;
    end


    % Data validation
    if size(DataIn,1)~=size(IDColumn,1)
        error('size(data_columns,1)~=size(code_column,1)!');
    end

    fprintf('%s\n','Restructuring data...  ');
    for iSubject=1:size(DataIn,1)
        xASL_TrackProgress(iSubject,size(DataIn,1));
        % Check in which population the subject is
        index   = find( codes== IDColumn(iSubject) );
        if  ~isempty(index)  % skips NaNs
            iCode   = codes( index);
            if ~isnan(iCode)
                % Put the row data of this subject in a new row in the population matrix
                DataOut{iCode}(NextN(iCode),:,:,:,:,:,:,:,:,:)    = DataIn(iSubject,:,:,:,:,:,:,:,:,:);
                NextN(iCode)                                      = NextN(iCode)+1;
            end
        end
    end
    fprintf('%s\n',['matrix with ' num2str(size(DataIn,1)) ' subjects and X ' num2str(size(DataIn,2)) ' Y ' num2str(size(DataIn,3)) ' Z ' num2str(size(DataIn,4)) ' T ' num2str(size(DataIn,5)) ' was restructured into ' num2str(nCodes) ' datasets']);

end
