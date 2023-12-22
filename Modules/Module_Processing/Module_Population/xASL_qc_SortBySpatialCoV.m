function xASL_qc_SortBySpatialCoV(x, Threshold1, Threshold2)
%xASL_wrp_GetROIstatistics Compute statistics for each ROI
%
% FORMAT: xASL_qc_SortBySpatialCoV(x, Threshold1, Threshold2)
% 
% INPUT:
%   x                            - struct containing statistical pipeline environment parameters (REQUIRED)
%   Threshold1 & Threshold2      - spatial CoV thresholds that are used to bin the ASL images per:
%                                  CBF contrast < Threshold1 < vascular contrast < Threshold2 < artifactual contrast
%                                  (OPTIONAL, DEFAULT -> Threshold1 = 0.67; Thresholds2 = 1)
%
% OUTPUT: n/a
% -------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function organizes the ASL QC images in //analysis/Population/ASLCheck
%              into CBF, vascular and artifactual contrast per the
%              spatial CoV thresholds defined above, in folders:
%              //analysis/Population/ASLCheck/1_CBFContrast
%              //analysis/Population/ASLCheck/2_VascularContrast
%              //analysis/Population/ASLCheck/3_ArtifactualContrast
% 			   Invalid spatial CoV numbers (e.g. NaN) will go to:
% 			   //analysis/Population/ASLCheck/4_Unknown_sCoV
%              Note: these outputs need to be visually checked; but
%              pre-sorting them by spatial CoV puts them already in
%              categories that are quick to skim through and manually
%              correct by moving the images
%
%              The idea is then that only category 1 images are used
%              for perfusion (CBF) analyses, both categories 1 & 2 for the
%              vascular (sCoV) analyses, and the category 3 images excluded
%              from analysis.
%
%              PM: this code does not include multiple sessions per subject yet!
% 			   NB: this code uses the //analysis/Population/Stats/CoV_qCBF*TotalGM*.csv file,
% 				   make sure that this file isn't edited! 
% -------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_qc_SortBySpatialCoV(x);
% __________________________________
% Copyright 2015-2020 ExploreASL


%% Admin
if nargin<2 || isempty(Threshold1)
    Threshold1 = 0.67;
end
if nargin<3 || isempty(Threshold2)
    Threshold2 = 1;
end

%% Find spatial CoV stats
FileList = xASL_adm_GetFileList(x.S.StatsDir, ['(?i).*CoV_qCBF.*TotalGM_n=' num2str(x.dataset.nSubjects) '_.*PVC0\.tsv$'], 'List',[0 Inf]);

if isempty(FileList)
    warning('Couldnt find spatial CoV information! File missing, skipping...');
    return;
end

% Now we take the most recent version
PathTSV = fullfile(x.S.StatsDir, FileList{end});

[~, CellTSV] = xASL_bids_csv2tsvReadWrite(PathTSV);
subjectList = CellTSV(3:end,1);

% If the subjectlist doesn't have sessions mentioned, we include them
if isempty(regexpi(subjectList{1}, 'ASL_\d+'))
    indexSession = find(strcmpi(CellTSV(1,:),'session'));
    if ~isempty(indexSession)
        sessionList = CellTSV(3:end,indexSession);
        subjectSessionList = strcat(subjectList,'_',sessionList);
    else
        subjectSessionList = subjectList;
        warning('Could not concatenate subjectlist and session list, this might go wrong');
    end
else
    subjectSessionList = subjectList;
end

sCoVList = CellTSV(3:end,end-2);

DirCBF = fullfile(x.D.ASLCheckDir, '1_CBFContrast');
DirVasc = fullfile(x.D.ASLCheckDir, '2_VascularContrast');
DirArti = fullfile(x.D.ASLCheckDir, '3_ArtifactContrast');
DirUnknown_sCoV = fullfile(x.D.ASLCheckDir, '4_Unknown_sCoV');

xASL_delete(DirCBF, true);
xASL_delete(DirVasc, true);
xASL_delete(DirArti, true);
xASL_delete(DirUnknown_sCoV, true);

%% Move the images
fprintf('Sorting ASLCheck QC images for spatial CoV:   ');
for iSubject=1:x.dataset.nSubjects
    for iSession=1:x.dataset.nSessions
        iSubjSess = (iSubject-1)*x.dataset.nSessions+iSession;
        NameSubjSess = [x.SUBJECTS{iSubject} '_' x.SESSIONS{iSession}];
        xASL_TrackProgress(iSubjSess, x.dataset.nSubjects * x.dataset.nSessions);
    
        % find current JPG
        JPGList = xASL_adm_GetFileList(x.D.ASLCheckDir, ['(?i)^Tra_qCBF_' NameSubjSess '.*'],'List', [0 Inf]);
        % find current sCoV
        Index = find(cellfun(@(y) strcmp(y, NameSubjSess), subjectSessionList));
        
        if isempty(JPGList) % if we did not find an ASL jpg image, skip this subject-session
            fprintf('\n');
            warning(['Missing transversal ASL jpg image in ' x.D.ASLCheckDir]);
            fprintf('%s\n\n', ['Did something go wrong in the ASL module for ' NameSubjSess '?']);
        else
            if isempty(Index) % if we cannot determine the sCoV, copy this subject-session to 4_Unknown_sCoV
                fprintf('\n');
                warning(['Missing sCoV value in ' PathTSV]);
                fprintf('%s\n\n', ['Did something go wrong in the ASL module for ' NameSubjSess '?']);
                Ddir = DirUnknown_sCoV;
            else
                sCoV = xASL_str2num(sCoVList{Index});

                % find folder
                if ~isfinite(sCoV) || ~isnumeric(sCoV) || sCoV<=0
                    % for illegal sCoV, 4_Unknown_sCoV
                    Ddir = DirUnknown_sCoV;
                elseif sCoV<Threshold1
                    % for relatively low sCoV, 1_CBFContrast
                    Ddir = DirCBF;
                elseif sCoV<Threshold2
                    % for intermediate sCoV, 2_VascularContrast
                    Ddir = DirVasc;
                else % for high sCoV, 3_ArtifactContrast
                    Ddir = DirArti;
                end
            end
    
            % copy files
            for iJPG=1:length(JPGList)
                OriFile = fullfile(x.D.ASLCheckDir, JPGList{iJPG});
                DestFile = fullfile(Ddir, JPGList{iJPG});
                xASL_adm_CreateDir(Ddir);
                xASL_Copy(OriFile, DestFile, true);
            end
        end
    end
end
fprintf('\n');

end
    