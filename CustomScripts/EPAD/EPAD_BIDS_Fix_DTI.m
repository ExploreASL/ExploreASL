function EPAD_BIDS_Fix_DTI(AnalysisDir)
%EPAD_BIDS_Fix_DTI Fix dcm2nii DTI conversion per BIDS
%
% FORMAT: EPAD_BIDS_Fix_DTI(AnalysisDir)
%
% INPUT:
%   AnalysisDir - path to folder containing the NIfTI/BIDS data (e.g. /data/RAD/share/EPAD/analysis)
%
% OUTPUT: n/a
% OUTPUT FILES:
% *NORMPE.nii
% including json sidecars
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function manages the DTI BIDS structure, for all subjects.
%              1) split the ADC & DTI source data (some vendors have this concatenated in a single NIfTI)
%              2) create the NormPE NIfTI from the first B0 volume from the source data
%                 the NormPE (i.e. the normal phase encoding image) is needed for TopUp distortion correction
%                 this file is not scanned separately as the DTI sequence
%                 usually also contains a B0 image. This function assumes
%                 that there is a .bval file with the same folder &
%                 filename as the dwi.nii, containing the b values. It will
%                 search for the first B0 image & copy this (including the
%                 _parms.mat & json sidecars) as a separate _NormPE.nii
%                 file.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: EPAD_BIDS_Fix_DTI(AnalysisDir);
% __________________________________
% Copyright 2015-2019 ExploreASL

SubjectList = xASL_adm_GetFsList(AnalysisDir,'^\d{3}EPAD\d*(|_\d*)$', true, [], [], [0 Inf]);

if isempty(SubjectList)
    warning('Didnt find subjects for DTI curation, skipping');
    return;
end

fprintf('%s','Adjusting DTI NIfTIs:  0%');

for iSubject=1:length(SubjectList)
    xASL_TrackProgress(iSubject,length(SubjectList));
    % find the first dwi file
    ScanDir = fullfile(AnalysisDir, SubjectList{iSubject},'dwi_1');
    Filelist = sort(xASL_adm_GetFileList(ScanDir, '^dwi(|_dwi)(|_run(-|_)\d*).*\.(nii|nii\.gz)$', 'FPList', [0 Inf]));

    %% FIRST WE CHECK THE SEPARATE ADC IMAGE (should go in the Reconstructions Dir?)
    if length(Filelist)>1 && ~isempty(strfind(lower(Filelist{2}),'adc'))
        tNii1 = xASL_io_Nifti2Im(Filelist{1});
        tNii2 = xASL_io_Nifti2Im(Filelist{2});

        if size(tNii2,4)==size(tNii1,4)+1 % assuming the last frame is the ADC frame
            [Fpath, Ffile, Fext] = xASL_fileparts(Filelist{1});
            ADCPath = fullfile(Fpath, [Ffile '_ADC' Fext]);
            xASL_io_SaveNifti(Filelist{1}, ADCPath, tNii2(:,:,:,end), [], 0);
            xASL_delete(Filelist{2});

            OriJSON = fullfile(Fpath, [Ffile '.json']);
            DestJSON = fullfile(Fpath, [Ffile '_ADC.json']);
            xASL_Copy(OriJSON, DestJSON);
        end
    end



    %% THEN WE CREATE THE NORMPE NIFTI
    if ~isempty(Filelist)
        PathNii = Filelist{1};
        [Fpath, Ffile, Fext] = xASL_fileparts(PathNii);
        PathBVAL = fullfile(Fpath, [Ffile '.bval']); % find the .bval file

        if exist(PathBVAL,'file')
            L = load(PathBVAL);
            IndexB0 = min(find(L==0));
        else
            warning(['Could not find ' PathBVAL, ': check contrast NormPE']);
            IndexB0 = 1;
        end

        if xASL_exist(PathNii,'file')
            % save NormPE NIfTI
            tNii = xASL_io_Nifti2Im(PathNii);
            ImNormPE = tNii(:,:,:,IndexB0);
            SavePathNii = fullfile(Fpath, [Ffile '_NormPE' Fext]);
            
            try
                xASL_io_SaveNifti(PathNii, SavePathNii, ImNormPE, [], 0);
            catch
                warning(['Broken NIfTI? check: ' PathNii]);
                continue;
            end

            % Copy JSON, issue warning if not existing
            PathJSON = fullfile(Fpath, [Ffile '.json']);
            SavePathJSON = fullfile(Fpath, [Ffile '_NormPE.json']);
            if exist(PathJSON, 'file')
                xASL_Copy(PathJSON, SavePathJSON, true);
            else
                warning(['Could not find ' PathJSON]);
            end
            
            % Copy _parms.mat if exist (but this is phasing out, don't
            % issue warning here)
            
            PathMAT = fullfile(Fpath, [Ffile '_parms.mat']);
            SavePathMAT = fullfile(Fpath, [Ffile '_NormPE_parms.mat']);
            if exist(PathMAT,'file')
                xASL_Copy(PathMAT, SavePathMAT, true);
            else
                PathMAT = xASL_adm_GetFileList(Fpath, 'dwi.*_parms\.mat','FPList',[0 Inf]);
                if ~isempty(PathMAT)
                    xASL_Copy(PathMAT{1}, SavePathMAT, true);
                end
            end
        end
    end
end

xASL_TrackProgress(1, 1);
fprintf('\n');

end
