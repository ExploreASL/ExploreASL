function regroupScripts(OriDir, DstDir)
% Regroups the data from subject/sequence/asl structure to 
% sequence/subject/asl
%
% INPUT:
%     OriDir - Original directory where the data is stored (REQUIRED)
%     DstDir - Destination directory for the newly regrouped data (OPTIONAL, DEFAULT = OriDir)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Regroups data files to comply with further standard
%               processing. The original data structure is of supposed to
%               be of the following subject/sequence/asl. Converts the
%               structure to sequence/subject/asl
%
% EXAMPLE: regroupScripts(/mnt/c/Backup/ASL/StudyName/analysis, /mnt/c/Backup/ASL/StudyName/analysis_2, '.*_2.*, 'true); Move all second time point folders & files to a new location, keeping directory structure intact
%
% 2020-12-13 Yeva Prysiazhniuk. This file was created as a part of ExploreASL
% See LICENSE for details.

%OriDir = 'C:\Users\yevap\Documents\CVUT\Diplomka-Petr\pet\projekte\asl\originalDataZipped\Philips\data';
%DstDir = 'C:\Users\yevap\Documents\CVUT\Diplomka-Petr\pet\projekte\asl\originalDataZipped\Philips\data_regrouped';

if nargin<2
    % Destinationa directory is created in the original directory - to be
    % changed in the future
    mkdir(OriDir, 'regrouped');
    DstDir = strcat(OriDir, '\regrouped');
end

files = dir(OriDir);
subjects = files(~ismember({files.name},{'.','..','PaxHeader'}) & [files.isdir]);
subDirs = [strcat({subjects.folder},'\',{subjects.name})];

seqNames = {};
for j = 1:length(subDirs)
    currSubName = subjects(j).name;
    
    seqDirFrom = dir(strcat(subDirs{j},'\','analysis'));
    seqDirFrom = seqDirFrom(~ismember({seqDirFrom.name},{'.','..','PaxHeader'}) & [seqDirFrom.isdir]);
    
    subT1From = strcat(seqDirFrom(1).folder, '\', seqDirFrom(1).name, '\T1w_1\');
    
    for k = 1:numel(seqDirFrom)
        currSeqName = seqDirFrom(k).name;
        
        if ~any(strcmp(seqNames,currSeqName))
            % If the folder is not created yet
            seqNames{end+1} = currSeqName;
            mkdir(DstDir, currSeqName);
        end
        
        subDirTo = dir(strcat(DstDir, '\', currSeqName));
        mkdir(subDirTo(1).folder, currSubName);
        dataDirTo = strcat(subDirTo(j).folder, '\', currSubName);
        dataDirFrom = strcat(seqDirFrom(k).folder, '\', currSeqName);
        
        % Copy T1s from first sequence
        copyfile(strcat(subT1From, 'T1w.*'), dataDirTo);
        
        % Copy ASL in a separate folder
        copyfile(strcat(dataDirFrom, '\ASL_1'), strcat(dataDirTo, '\ASL_1'));
        
    end
    
end