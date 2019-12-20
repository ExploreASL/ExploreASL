function [TSVfile, TempCell] = xASL_adm_csv2tsv(InPath, IfDelete, IfCreateTSV)
%xASL_adm_csv2tsv Convert a csv to a tsv file, for BIDS compliance.
% Outputs also as cell.
% Can also take tsv as input, in which case it only passes out to cell

% InPath = 'C:\Backup\ASL\Sleep_2018\analysis\dartel\StatsMaps\Sleep_2018_Test\permute_LongitudinalTimePoint\ClusterStats_Sleep_2018_Test 1w ws-ANOVA LongitudinalTimePoint.csv';

if ~exist('IfDelete', 'var')
    IfDelete = true;
end

if ~exist('IfCreateTSV', 'var')
    IfCreateTSV = true;
end

[Fpath, Ffile, Fext]     = fileparts(InPath);

if ~exist(InPath, 'file') && strcmp(Fext, '.csv')
    % first try replacing the extension by csv, as this script converts everything to tsv
    Fext = '.tsv';
    TSVPath = fullfile(Fpath, [Ffile '.tsv']);
    if exist(TSVPath, 'file')
        InPath = TSVPath;
    else
        error(['Could not find ' InPath]);
    end
end

if      strcmp(Fext,'.tsv')
        
        ReadCell 	= textread(InPath, '%s', 'delimiter','\t\n', 'bufsize',10000000);
        % now we have a 1D cell structure, where the \n became empty cells (other stuff is skipped)
        iRow    = 1;
        iColumn = 1;
        for iT=1:length(ReadCell)
            if  isempty(ReadCell{iT})
                iRow        = iRow+1;
                iColumn     = 1;
            else
                TempCell{iRow,iColumn}  = ReadCell{iT};
                iColumn     = iColumn+1;
            end
        end
            
   
elseif  strcmp(Fext,'.csv')

    TempText 	= textread(InPath, '%s', 'delimiter','\n', 'bufsize',10000000);

    for iT=1:size(TempText,1)
        CommaIndices = find(TempText{iT,1}==',');
        CommaIndices = [0 CommaIndices length(TempText{iT,1})+1];
        for iC=1:length(CommaIndices)-1
            CI1                 = CommaIndices(iC)+1;
            CI2                 = CommaIndices(iC+1)-1;
            TempCell{iT,iC}     = TempText{iT,1}(CI1:CI2);

            %% Replaces spaces by underscores
            iSpace              = find(TempCell{iT,iC}==' ');
            while ~isempty(iSpace)
                if      iSpace(1)==1
                        TempCell{iT,iC} = TempCell{iT,iC}(2:end);
                elseif  iSpace(1)==length(TempCell{iT,iC})
                        TempCell{iT,iC} = TempCell{iT,iC}(1:end-1);
                else
                    TempCell{iT,iC} = [TempCell{iT,iC}(1:iSpace(1)-1) '_' TempCell{iT,iC}(iSpace(1)+1:end)];
                end
                iSpace              = find(TempCell{iT,iC}==' ');
            end

        end
    end
end

% [numericData, textData, rawData] = xlsread(InPath); doesnt work on all platforms

% TempTable = readtable(TempFile,'Delimiter',',','format','%s','ReadVariableNames',0); 
% use detectImportOptions here to get correct format
% if we need to avoid xlsread (not available on all platforms?)

TSVfile     = fullfile(Fpath,[Ffile '.tsv']);

if ~exist(TSVfile,'file') && IfCreateTSV

    FID         = fopen(TSVfile ,'w');

    for iX=1:size(TempCell,1)
        for iY=1:size(TempCell,2)
            fprintf(FID,'%s\t',xASL_num2str(TempCell{iX,iY}));
        end
        fprintf(FID,'\n');
    end

    fclose(FID);
end

if  IfDelete && exist(TSVfile,'file') && strcmp(InPath(end-2:end),'csv')
    xASL_delete(InPath);
end
    
end

