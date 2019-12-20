function xASL_io_SplitASL_M0(InPath,iM0)
%xASL_io_SplitASL_M0 Splits ASL & M0 images, when they are within the same NIfTI
%
% FORMAT: xASL_io_SplitASL_M0(InPath,iM0)
%
% INPUT:
%   InPath      - path to ASL NIfTI file (e.g. //analysis/ASL_1/ASL4D.nii) (REQUIRED)
%   iM0         - index/indices of volume(s) containing the M0 (REQUIRED)
%
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function splits ASL4D & M0 if they were in the same sequence.
%              If dcm2niiX has already splitted the ASL4D NIfTI, it is
%              reconstructed first.
%
%              Vendor product sequence examples:
%              GE 3D spiral sometimes puts the M0 at the last volume of the series -> iM0 = [2];
%              Philips 3D GRASE puts the M0 as control-label volume pair -> iM0 = [1 2];
%              Siemens 3D GRASE puts the M0 as the first volume -> iM0 = 1;
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:
% xASL_io_SplitASL_M0('/data/RAD/share/EPAD500/010EPAD00001/ASL_1/ASL4D.nii', [1 2]);
% __________________________________
% Copyright 2015-2019 ExploreASL


    [Fpath, Ffile]  = xASL_fileparts(InPath);

    % Split_ASL_M0
    ASLname = fullfile(Fpath, 'ASL4D.nii');
    BackupName = fullfile(Fpath,'ASL4D_Source.nii.gz');

    ASLlist = xASL_adm_GetFileList(Fpath, ['^' Ffile '(|_.)'  '(|_\d)' '\.nii$'], 'FPList', [0 Inf]);
  
    if isempty(ASLlist)
        error([ASLname ' didnt exist, skipping']);
    end
    
    if ~xASL_exist(BackupName) % otherwise was already split

        [Fpath, Ffile] = xASL_fileparts(ASLlist{1});

        OriMATPath = fullfile(Fpath,[Ffile '_parms.mat']);
        OriJSONPath = fullfile(Fpath,[Ffile '.json']);

        ASLMATPath = fullfile(Fpath, 'ASL4D_parms.mat');
        ASLJSONPath = fullfile(Fpath, 'ASL4D.json');
        M0MATPath = fullfile(Fpath, 'M0_parms.mat');
        M0JSONPath = fullfile(Fpath, 'M0.json');
        BackupMATPath = fullfile(Fpath, 'ASL4D_Source_parms.mat');
        BackupJSONPath = fullfile(Fpath, 'ASL4D_Source.json');

        %% First concatenate NIfTIs
        if length(ASLlist)>1
            % Reconstruct the ASL4D first
            for iNii=1:length(ASLlist)
                tnii = xASL_io_ReadNifti(ASLlist{iNii});
                RescaleSlope(iNii) = tnii.dat.scl_slope;
                
                tIM = xASL_io_Nifti2Im(tnii.dat(:,:,:,:,:,:,:));
                Dim4 = size(tIM,4);
                if iNii==1
                    TotIm = tIM;
                else
                    TotIm(:,:,:,CDim4+1:CDim4+Dim4) = tIM;
                end
                CDim4 = size(TotIm,4);
            end
            
            if length(unique(RescaleSlope))>1
                warning('Rescaleslopes were not the same between concatenated NIfTIs, skipping...');
                return;
            else
                % Save concatenated ASL series as backup ASL
                xASL_io_SaveNifti(ASLlist{1}, BackupName, TotIm, [], false);
                % 
                if exist(OriMATPath,'file') && ~strcmp(OriMATPath, BackupMATPath)
                    xASL_Move(OriMATPath, BackupMATPath);
                end
                if exist(OriJSONPath,'file') && ~strcmp(OriJSONPath, BackupJSONPath)
                    xASL_Move(OriJSONPath, BackupJSONPath);
                end
                
                for iNii=1:length(ASLlist)
                    xASL_delete(ASLlist{iNii});
                    [Fpath, Ffile] = xASL_fileparts(ASLlist{iNii});
                    JSON2Delete = fullfile(Fpath, [Ffile '.json']);
                    MAT2Delete = fullfile(Fpath, [Ffile '_parms.mat']);
                    xASL_delete(JSON2Delete);
                end
            end
            
        else % backup the ASL4D.nii & sidecars
            xASL_Move(ASLlist{1}, BackupName, true); 
            if exist(OriMATPath,'file') && ~strcmp(OriMATPath, BackupMATPath)
                xASL_Move(OriMATPath, BackupMATPath);
            end
            if exist(OriJSONPath,'file') && ~strcmp(OriJSONPath, BackupJSONPath)
                xASL_Move(OriJSONPath, BackupJSONPath);
            end
        end

        %% Save M0 NIfTI
        Path_M0 = fullfile(Fpath,'M0.nii');
        if ~xASL_exist(Path_M0,'file') % don't overwrite
            tIM = xASL_io_Nifti2Im(BackupName);
            xASL_io_SaveNifti(BackupName,Path_M0 ,tIM(:,:,:,iM0),[],0);
            bCreateM0 = true;
        else
            bCreateM0 = false;
        end

        %% Determine ASL indices
        ASLindices          = 1:1:size(tIM,4);
        IndexASL            = ones(1,size(tIM,4));
        IndexASL(iM0)       = 0;
        ASLindices          = ASLindices(logical(IndexASL));
        % Check for even number of volumes
        nASL                = numel(ASLindices);
        IsEven              = nASL/2==round(nASL/2); % even can be divided by 2
        ContainsTimeSeries  = nASL>5;
        if ContainsTimeSeries && ~IsEven
            warning('After removing M0 volume, no even number of volumes (control-label)');
        end

        %% Save ASL4D NIfTI
        xASL_io_SaveNifti(BackupName,ASLname,tIM(:,:,:,ASLindices),[],false);

        %% Copy sidecars
        if exist(BackupMATPath,'file') && bCreateM0
            xASL_Copy(BackupMATPath, M0MATPath);
        end
        if exist(BackupMATPath,'file') && ~strcmp(BackupMATPath, ASLMATPath)
            xASL_Copy(BackupMATPath, ASLMATPath);
        end
        if exist(BackupJSONPath,'file') && bCreateM0
            xASL_Copy(BackupJSONPath, M0JSONPath);
        end
        if exist(BackupJSONPath,'file') && ~strcmp(BackupJSONPath, ASLJSONPath)
            xASL_Copy(BackupJSONPath, ASLJSONPath);
        end
    end

end
