function [UnequalList] = UnitTest_xASL_io_DcmtkRead(DcmPath)
%UnitTest_xASL_io_DcmtkRead Unit test xASL_io_DcmtkRead

% Used for testing
% DcmPath = 'C:\Backup\ASL\EPAD\raw\040\040EPAD7230';

DcmList = xASL_adm_GetFileList(DcmPath, '^.*[^.(csv|tsv|gz|json)]$','FPListRec',[0 Inf]); % \.dcm|.*\.img|.*\.IMA|
UnequalList = struct;
Ind1 = length(DcmPath)+1;

fprintf('Performing unit test on DcmtkRead:   ');

for iD=1:length(DcmList)
    xASL_TrackProgress(iD,length(DcmList));
    TKRead = xASL_io_DcmtkRead(DcmList{iD});
    Info = dicominfo(DcmList{iD});
    
    DICOMfields = fields(TKRead);
    
    CurrFieldN = num2str(length(fields(UnequalList)));
    
    for iL=1:length(DICOMfields)
        if ~isfield(Info,DICOMfields{iL}) && isempty(TKRead.(DICOMfields{iL}))
            % this field didnt exist in the DICOM file, same answer for both dicominfo & dcmtk
        elseif isfield(Info,DICOMfields{iL}) && ~isempty(TKRead.(DICOMfields{iL}))
            % both codes read this field
            
            ContentTK = TKRead.(DICOMfields{iL});
            ContentInfo = Info.(DICOMfields{iL});          
            
            if ~isequal(ContentTK, ContentInfo)
                NewField = ['Unequal_field_' DICOMfields{iL}];
                if ~isfield(UnequalList,NewField)
                    IndH = 1;
                else
                    IndH = length(UnequalList.(NewField));
                end
                UnequalList.(NewField){IndH,1} = ['TK: ' xASL_num2str(ContentTK) ', dicominfo: ' xASL_num2str(ContentTK)];
                UnequalList.(NewField){IndH,2} = DcmList{iD}(Ind1:end);
            end

        elseif ~isfield(Info,DICOMfields{iL}) && ~isempty(TKRead.(DICOMfields{iL}))
                NewField = ['Missing_field_' DICOMfields{iL}];
                if ~isfield(UnequalList,NewField)
                    IndH = 1;
                else
                    IndH = length(UnequalList.(NewField));
                end
                UnequalList.(NewField){IndH,1} = 'Missing from Dicominfo';
                UnequalList.(NewField){IndH,2} = DcmList{iD}(Ind1:end);
            
        elseif isfield(Info,DICOMfields{iL}) && isempty(TKRead.(DICOMfields{iL}))
                NewField = ['Missing_field_' DICOMfields{iL}];
                if ~isfield(UnequalList,NewField)
                    IndH = 1;
                else
                    IndH = length(UnequalList.(NewField));
                end
                UnequalList.(NewField){IndH,1} = 'Missing from Dicominfo';
                UnequalList.(NewField){IndH,2} = DcmList{iD}(Ind1:end);            
            
            
            
            UnequalList.(['Missing_field_' DICOMfields{iL}]) = 'Missing from TKRead';
            UnequalList.(['DICOMpath_' CurrFieldN]) = DcmList{iD}(Ind1:end);
        end
    end
   
    % Save the list
    SavePath = fullfile(DcmPath,'UnequalList.json');
    spm_jsonwrite(SavePath, UnequalList);
end

fprintf('\n');

end

% TODO: script needs to check whether it is a number or string, for a number, round it to find equal content