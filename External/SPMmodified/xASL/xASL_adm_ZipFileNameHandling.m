function [srcOut, dstOut] = xASL_adm_ZipFileNameHandling(srcIn, dstIn, bDeleteOldest)
% xASL_adm_ZipFileNameHandling Works only for NIFTI files - changes the source and dest filenames according to what are the available inputs
%
% FORMAT: [srcOut, dstOut] = xASL_adm_ZipFileNameHandling(srcIn, dstIn)
%
% INPUT:
%   srcIn         - Source file on input (REQUIRED)
%   dstIn         - Destination file on input (REQUIRED)
%   bDeleteOldest - boolean to specify whether oldest file should be
%                   deleted if both .nii & .nii.gz counterparts exist, but
%                   are unequal (OPTIONAL, DEFAULT=false). WARNING PUT THIS
%                   ON TRUE AT OWN RISK
% OUTPUT:
%   srcOut   - Source file on output
%   dstOut   - Destination file on output
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Adjusts the source and destination filenames of a nifti file to reflect if NII or NII.GZ exist on the input.
%              If either .nii or .nii.gz is corrupt, it automatically deletes the corrupt one and keeps the healthy one, 
%              while reporting a warning. This happens when you restart the pipeline after it crashed, if it crashed while unzipping. 
%
% EXAMPLE: xASL_adm_ZipFileNameHandling('c:\User\path\file.nii.gz', 'c:\User\path2\file.nii');
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright (C) 2015-2019 ExploreASL
%
% 2018-01-01 JP, HM

    % Admin
    if nargin<2 || isempty(dstIn)
        dstIn = srcIn;
    end
    if nargin<3 || isempty(bDeleteOldest)
        bDeleteOldest = false;
    end

    % Break to path/file/extension
    [srcPath, srcFile, srcExt] = xASL_fileparts(srcIn);
    [dstPath, dstFile] = xASL_fileparts(dstIn);

	% If nothing changes, keep the output same as the input
	dstOut = dstIn;
	srcOut = srcIn;

	% Checks if the source concerns a .nii or .nii.gz
    if ~isempty(findstr(srcExt,'.nii'))

        srcFileNII = fullfile(srcPath,[srcFile '.nii']);
        srcFileGZ = fullfile(srcPath,[srcFile '.nii.gz']);

        % If the temporary file also exists, delete it first
        if exist([srcFileGZ '.tmp'], 'file')
            delete([srcFileGZ '.tmp']);
        end
        
		% If both input files exist, then try to delete one of them
		if exist(srcFileNII,'file') && exist(srcFileGZ,'file')
			% Unzip the GZ file to a temporary directory
			gunzip(srcFileGZ, [srcFileGZ '.tmp']);

            % Load both the original and new unzipped
            % If one is corrupt, replace it by the other
            GZcounterpart = fullfile([srcFileGZ '.tmp'],[srcFile '.nii']);
            try % first we check if the .nii is corrupt
                srcDatOrig = nifti(srcFileNII);
                srcImOrig = xASL_io_Nifti2Im(srcDatOrig);
            catch ME
                % This file is corrupt, try to fix with the other
                warning(ME.message);
                fprintf('%s\n', ['Corrupt NIfTI: ' srcFileNII]);
                fprintf('%s\n', 'Trying to replace by GZ counterpart');
                clear srcDatOrig srcImOrig
                xASL_SysCopy(GZcounterpart, srcFileNII, true);
                srcDatOrig = nifti(srcFileNII);
                srcImOrig = xASL_io_Nifti2Im(srcDatOrig);
            end
            try % secondly we check if the .nii.gz is corrupt
                srcDatZipd = nifti(GZcounterpart);
                srcImZipd = xASL_io_Nifti2Im(srcDatZipd);
            catch ME
                % This file is corrupt, try to fix with the other
                warning(ME.message);
                fprintf('%s\n', ['Corrupt NIfTI: ' GZcounterpart]);
                clear srcDatZipd srcImZipd
                xASL_SysCopy(srcFileNII, GZcounterpart, true);
                srcDatZipd = nifti(GZcounterpart);
                srcImZipd = xASL_io_Nifti2Im(srcDatZipd);
            end

			% Delete the temporary file
			delete(fullfile([srcFileGZ '.tmp'],'*'));
			rmdir([srcFileGZ '.tmp']);

			% Compare files:
            % First we compare the NIfTI header inside the hrd subfield
			AreEqual = true;
            DifferIn = '';
            FieldList = fields(srcDatOrig.hdr);
			for iL=1:length(FieldList)
				if ~isfield(srcDatZipd.hdr,(FieldList{iL})) || ~isequal(srcDatOrig.hdr.(FieldList{iL}),srcDatZipd.hdr.(FieldList{iL}))
					AreEqual = false;
					if ~isempty(DifferIn); DifferIn = [DifferIn '_']; end
					DifferIn = [DifferIn '.hdr.' FieldList{iL}];
				end
            end
            if ~isempty(DifferIn)
                fprintf('.nii & .nii.gz counterparts differed in:\n');
                fprintf('%s\n', xASL_num2str(DifferIn));
            end
            
            % Second, we compare the image matrices
			if ~isequal(srcImOrig,srcImZipd)
				AreEqual = false;
				if ~isempty(DifferIn); DifferIn = [DifferIn '_']; end
				DifferIn = [DifferIn 'image'];
                fprintf('.nii & .nii.gz counterparts differed in their image matrix\n');
			end

           % Third, we compare the orientation matrices (header outside the hdr subfield)
			FieldList = {'mat';'mat0';'mat_intent';'mat0_intent';'timing';'descrip'};

			for iL=1:length(FieldList)
				if isfield(srcDatOrig,FieldList{iL})
					if ~isfield(srcDatZipd,FieldList{iL}) || ~isequal(srcDatOrig.(FieldList{iL}),srcDatZipd.(FieldList{iL}))
						AreEqual = false;
						if ~isempty(DifferIn); DifferIn = [DifferIn '_']; end
						DifferIn = [DifferIn '.' FieldList{iL}];
					end
				end
			end

			if AreEqual
				% If same - delete one
				if strcmp(srcFileGZ,srcIn)
					% Preferentially keep the one specified in the input
					delete(srcFileNII);
				else
					delete(srcFileGZ);
				end
            elseif bDeleteOldest % if unequal
                FileInfoNII = dir(srcFileNII);
                FileInfoGZ = dir(srcFileGZ);
                warning('Found two unequal valid NIfTI(.gz) files with same name, deleting oldest one');
                % Delete the oldest one
                if FileInfoNII.datenum>FileInfoGZ.datenum
                    delete(srcFileGZ);
                    fprintf('Deleted .GZ counterpart\n');
                else
                    delete(srcFileNII);
                    fprintf('Deleted .NII counterpart\n');
                end
            else
				% Otherwise, throw an error
				error(['Two files with same name, but different ' DifferIn ': \n %s \n %s'], srcFileNII, srcFileGZ);
			end
		end

        % Now only if .nii or .nii.gz doesnt exist, use the existing one
        if ~exist(srcFileNII,'file') && exist(srcFileGZ,'file')
            srcOut = srcFileGZ;
            srcExt = '.nii.gz';
        elseif ~exist(srcFileGZ,'file') && exist(srcFileNII,'file')
            srcOut = srcFileNII;
            srcExt = '.nii';
        end

        % Here, to keep it simple, we force DestExt to be the same as SrcExt
        % This avoids unzipping/zipping, which doesnt have to be dealth with when copying
        dstOut = fullfile(dstPath,[dstFile, srcExt]);
    end

end
