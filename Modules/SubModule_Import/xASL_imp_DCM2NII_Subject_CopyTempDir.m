function xASL_imp_DCM2NII_Subject_CopyTempDir(nii_files, bClone2Source)
%xASL_imp_DCM2NII_Subject_CopyTempDir Make a copy of analysisdir in sourcedir.
%
% FORMAT: xASL_imp_DCM2NII_Subject_CopyTempDir(nii_files, bClone2Source)
% 
% INPUT:
%   nii_files     - List of NIfTI files (CELL ARRAY, REQUIRED)
%   bClone2Source - Clone to source (BOOLEAN, REQUIRED)
%
% OUTPUT:
%   n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Make a copy of temp dir in source dir.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     n/a
%
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Make a copy of temp dir in source dir
    if bClone2Source
        if ~isempty(nii_files)
            for iFile=1:length(nii_files)
                % replace 'temp' by 'source'
                [iStart, iEnd] = regexp(nii_files{iFile}, 'temp');
                DestPath = [nii_files{iFile}(1:iStart-1) 'source' nii_files{iFile}(iEnd+1:end)];
                xASL_Copy(nii_files{iFile}, DestPath, true);
                % do the same for other extensions
                Extensions = {'.json' '_parms.json'};
                for iExt=1:length(Extensions)
                    [Fpath, Ffile] = xASL_fileparts(nii_files{iFile});
                    CopyPath = fullfile(Fpath, [Ffile Extensions{iExt}]);
                    [Fpath, Ffile] = xASL_fileparts(DestPath);
                    DestPath = fullfile(Fpath, [Ffile Extensions{iExt}]);
                    if xASL_exist(CopyPath)
                        xASL_Copy(CopyPath, DestPath, true);
                    end
                end
            end
        end
    end

end

