function [niifiles, ScanNameOut, usedinput, msg] = xASL_io_dcm2nii(inpath, destdir, series_name, varargin)
% xASL_io_dcm2nii
% Convert DICOM NIfTI/BIDS format using the dcm2nii command line utility
% (http://www.nitrc.org/projects/mricron)
%
% - inpath           path to dicom folder, dicom file, PAR-file or REC-file. In case of a dicom folder, 'DcmExt' will be used to
%                    filter the folder contents. In case of a dicom file, all dicom files in the same folder will be used.
% - 'Verbose'        optional: set to true or false to switch terminal feedback; default false
% - 'Overwrite'      optional: set to true or false to overwrite existing files; default false
% - 'DicomFilter'    optional: regular expression used to find dicom files; default '^.+\.dcm$'
% - 'Keep'           optional: vector with indices of output files to keep; default [1:Inf] (all files)
% - 'IniPath'        optional: string with path of dcm2nii configuration file; default is to use ./mricron/dcm2nii-custom.ini
% - 'ExePath'        optional: string with path of dcm2nii application; default is to use internal copy in ./mricron folder
% - 'Version'        optional: string with version stamp of dcm2nii application (yyyymmdd: 20101105, 20130606); default is 20101105
%
% The first DICOM file will be used if inpath is a directory.
%   Detailed explanation goes here

    niifiles = {};
    msg = [];

    p = inputParser;
    addOptional(p, 'Verbose', false, @islogical);
    addOptional(p, 'Keep', [], @isvector);
    addOptional(p, 'Overwrite', false, @islogical);
    addOptional(p, 'DicomFilter', '^.+\.dcm$', @ischar);
    addOptional(p, 'ExePath', [], @ischar);
    addOptional(p, 'IniPath', [], @ischar);
    addOptional(p, 'Version', '20101105', @ischar);
    addOptional(p, 'x', struct, @isstruct);

    %% parse parameters
    parse(p,varargin{:});
    parms = p.Results;

    %% locate dcm2nii executable
    if ismac()
        parms.Version = '20181125'; % mac is incompatible with older versions
    end

    if ~isfield(parms,'x')
        warning('Cannot retrieve ExploreASL folder from x-structure, using current path instead');
        MyPath = pwd;
    else
        MyPath = parms.x.MyPath;
    end

    mricron_path = fullfile(MyPath,'External','MRIcron');
    if isempty(parms.ExePath)
        mricron_version_path = fullfile(mricron_path, parms.Version);
        switch computer()
            case {'PCWIN', 'PCWIN64'}
                parms.ExePath = fullfile(mricron_version_path,'dcm2nii.exe');
            case 'GLNX86'
                parms.ExePath = fullfile(mricron_version_path,'dcm2nii-lx32');
            case 'GLNXA64'
                parms.ExePath = fullfile(mricron_version_path,'dcm2nii-lx64');
            case {'MACI','MACI64'}
                parms.ExePath = fullfile(mricron_version_path,'dcm2nii-osx');
            otherwise
                error('Unknown computer architecture: %s',computer());
        end
    end
    [stat, attr] = fileattrib(parms.ExePath);
    if ~stat
        error('dcm2nii application not found: %s',parms.ExePath);
    end
    if ~attr.UserExecute
        % try fixing
        try
            system(['chmod 777 ' parms.ExePath]);
        catch
            disp(attr)
            error('dcm2nii application is not executable!');
        end
    end

    %% set default arguments dcm2nii
%    '-a y -d n -e n -f y -g n -n y -p n -r n -v y -x n'; % -o will be appended below
    % -a Anonymize [remove identifying information]: Y,N = Y
    % -b load settings from specified inifile, e.g. '-b C:\set\t1.ini'
    % -c Collapse input folders: Y,N = N
    % -d Date in filename [filename.dcm -> 20061230122032.nii]: Y,N = Y
    % -e events (series/acq) in filename [filename.dcm -> s002a003.nii]: Y,N = Y
    % -f Source filename [e.g. filename.par -> filename.nii]: Y,N = Y
    % -g gzip output, filename.nii.gz [ignored if '-n n']: Y,N = Y
    % -i ID  in filename [filename.dcm -> johndoe.nii]: Y,N = Y
    % -n output .nii file [if no, create .hdr/.img pair]: Y,N = Y
    % -o Output Directory, e.g. 'C:\TEMP' (if unspecified, source directory is used)
    % -p Protocol in filename [filename.dcm -> TFE_T1.nii]: Y,N = Y
    % -r Reorient image to nearest orthogonal: Y,N
    % -s SPM2/Analyze not SPM5/NIfTI [ignored if '-n y']: Y,N = N
    % -v Convert every image in the directory: Y,N = Y
    % -x Reorient and crop 3D NIfTI images: Y,N = N

    %% set default arguments dcm2niiX
%   -1..-9 : gz compression level (1=fastest..9=smallest, default 6)
%   -b : BIDS sidecar (y/n/o [o=only: no NIfTI], default y)
%    -ba : anonymize BIDS (y/n, default y)
%   -c : comment stored in NIfTI aux_file (up to 24 characters)
%   -d : directory search depth. Convert DICOMs in sub-folders of in_folder? (0..9, default 5)
%   -f : filename (%a=antenna (coil) name, %b=basename, %c=comments, %d=description, %e=echo number, %f=folder name, %i=ID of patient, %j=seriesInstanceUID, %k=studyInstanceUID, %m=manufacturer, %n=name of patient, %p=protocol, %r=instance number, %s=series number, %t=time, %u=acquisition number, %v=vendor, %x=study ID; %z=sequence name; default '%f_%p_%t_%s')
%   -g : generate defaults file (y/n/o/i [o=only: reset and write defaults; i=ignore: reset defaults], default n)
%   -h : show help
%   -i : ignore derived, localizer and 2D images (y/n, default n)
%   -l : losslessly scale 16-bit integers to use dynamic range (y/n, default n)
%   -m : merge 2D slices from same series regardless of study time, echo, coil, orientation, etc. (y/n, default n)
%   -n : only convert this series number - can be used up to 16 times (default convert all)
%   -o : output directory (omit to save to input folder)
%   -p : Philips precise float (not display) scaling (y/n, default y)
%   -r : rename instead of convert DICOMs (y/n, default n)
%   -s : single file mode, do not convert other images in folder (y/n, default n)
%   -t : text notes includes private patient details (y/n, default n)
%   -v : verbose (n/y or 0/1/2 [no, yes, logorrheic], default 0)
%   -x : crop (y/n, default n)
%   -z : gz compress images (y/i/n/3, default n) [y=pigz, i=internal:miniz, n=no, 3=no,3D]

    %% Set dcm2niiX stuff
    if strcmp(parms.Version, '20181125')
        bBIDS = true;
        dcm2nii_args = ''; % -z y
    else
        bBIDS = false;
    end

    %% Check if we are reading a DICOM folder
    if exist(inpath,'dir')
        %% check if there are any DICOM files
        dicom_files = xASL_adm_GetFileList(inpath, parms.DicomFilter,'List');
        if isempty(dicom_files)
            error(['No dicom files match ' parms.DicomFilter ' in ' inpath ', skipping']);
            % this error will be catched in wrapper function, verbose & error catching will be dealt with there
        end

        %% sort the dicom filenames because it's nice to use the first (but not required!)
        dicom_files = sort(dicom_files);
        inpath = fullfile(inpath, dicom_files{1});
        bScanDicomFolder = true;
    else
        bScanDicomFolder = false;
    end

    if isempty(parms.IniPath)
        if bScanDicomFolder
            parms.IniPath = fullfile(mricron_path,'dcm2nii-dicom.ini');
        else
            parms.IniPath = fullfile(mricron_path,'dcm2nii-parrec.ini');
        end
    end

    if ~bBIDS
        dcm2nii_args = sprintf('-b "%s"', parms.IniPath);
    end


    %% check for existing targets
    bSingleExists = exist(fullfile(destdir, [series_name '.nii']),'file');
    bMultiExists = exist(fullfile(destdir, [series_name '_1.nii']),'file');
    if ~parms.Overwrite && (bSingleExists || bMultiExists)
        if bSingleExists
            niifiles{1} = fullfile(destdir, [series_name '.nii']);
        else
            niifiles = xASL_adm_GetFileList(destdir, ['^' regexptranslate('escape',series_name) '_\d+\.nii$'],'List');
            for iS=1:length(niifiles)
                niifiles{iS} = fullfile(destdir,niifiles{iS}); % must prefix path
            end
        end

        if nargout>1 % optionally return the used input file
            usedinput = inpath; % could be the original input or the dicom file used
        end
        ScanNameOut = series_name;

        if parms.Verbose; fprintf('%s\n', [ScanNameOut ' already existed, skipping...']); end
        return;
    end


    %% create temporary subfolder for converting
    temp_dir = fullfile(destdir, ['dcm2nii_temp_' series_name]);
%     if exist(temp_dir,'dir') && length(dir(temp_dir))>2
%         error('Temporary conversion folder exists and is not empty: %s', temp_dir);
%     end
%   IGNORE PREVIOUSLY STORED TEMPORARY FILES, THAT ARE STILL HERE BECAUSE OF EG CODE CRASHING
    if exist(temp_dir,'dir')
        xASL_adm_DeleteFileList(temp_dir,'.*', true, [0 Inf]);
    end
    xASL_adm_CreateDir(temp_dir);

    %% Run dcm2nii and move files to final destination using specified series name
    %% Use catch to remove temporary folder on error
    try
        %% execute
        if ispc()
            quote = '"';
        else
            quote = '''';
        end

        if ismac()
            cmdline = [parms.ExePath ' -o ' temp_dir ' ' inpath];
        else
            cmdline = [quote parms.ExePath quote ' ' dcm2nii_args ' -o ' quote temp_dir quote ' ' quote inpath quote];
        end

        if parms.Verbose
            fprintf('executing: [%s]\n',cmdline);
        end
        [status, msg] = system(cmdline);
        if parms.Verbose
            separatorline = repmat(char('-'),1,80);
            fprintf('%s\n%s%s\nstatus: %d\n',separatorline,msg,separatorline,status);
        end

        %% move/rename nifties to final destination
        niiEntries = xASL_adm_GetFileList(temp_dir, '.*\.nii$', 'FPListRec', [0 Inf]);
        if isempty(niiEntries)
            error('Empty output dcm2nii');
        else
            [niiPaths, niiNames, niiExt] = cellfun(@(y) xASL_fileparts(y), niiEntries, 'UniformOutput',false);
        end

        nNifties = length(niiNames);
        if nNifties==0
            error('Conversion failed: No nifti output files available');
        else
            if isempty(parms.Keep)
                parms.Keep = 1:nNifties;
            else
                parms.Keep = unique(parms.Keep);
                if max(parms.Keep)>nNifties
                    error('Invalid index: the requested volume (%d) exceeds the actual number of files (%d)',max(parms.Keep),nNifties);
                end
            end

            % Check if dcm2nii added some contrast information at the end of the FileName
            if length(parms.Keep)==1
                ContrastAdd{1} = [];
            else
                % first compare the filenames (for their equivalent length)
                EquiL = min(cellfun(@(y) length(y), niiNames));

                NN=1;
                for iN1=parms.Keep
                    for iN2=parms.Keep
                        CmpString(NN,:) = double(niiNames{iN1}(1:EquiL)~=niiNames{iN2}(1:EquiL));
                        NN=NN+1;
                    end
                end
                IndCommon = find(sum(CmpString,1)==size(CmpString,1)); % common filename differences (i.e. within equivalent filename length)
                for iN=parms.Keep
                    ContrastAdd{iN} = niiNames{iN}(IndCommon);
                    if length(niiNames{iN})>EquiL
                        ContrastAdd{iN} = [ContrastAdd{iN} niiNames{iN}(EquiL+1:end)];
                    end
                end
            end

            for iVolume=sort(parms.Keep)
                temp_file = niiEntries{iVolume};

                % Make BIDS compatible dest_file here
                DestFileName = series_name;

                [Ind1 Ind2] = regexp(niiNames{iVolume},'_run-\d*_','start','end');
                if ~isempty(Ind1) && ~isempty(Ind2) % add run_index suffix (if it is there)
                    RunName = niiNames{iVolume}(Ind1+1:Ind2-1);
                    DestFileName = [DestFileName '_' RunName];
                end
                if ~isempty(ContrastAdd{iVolume}) % add contrast suffix (if there is something that indicates it)
                    DestFileName = [DestFileName '_' ContrastAdd{iVolume}];
                end

                if iVolume==min(parms.Keep)
                    ScanNameOut = DestFileName;
                end

                if length(parms.Keep)>1 % add iVolume suffix (if there are multiple)
                    DestFileName = [DestFileName '_' int2str(iVolume)];
                end
                DestFileName = xASL_adm_CorrectName(DestFileName);

                dest_file = fullfile(destdir,[DestFileName '.nii']);
                xASL_Move(temp_file, dest_file, parms.Overwrite, parms.Verbose);
                niifiles{end+1} = dest_file; %#ok<AGROW>

                % Do the same for the BIDS files
                BIDSext = {'.json' '.bval' '.bvec'};
                for iB=1:length(BIDSext)
                    [Gpath, Gfile] = xASL_fileparts(temp_file);
                    temp_BIDS = fullfile(Gpath, [Gfile BIDSext{iB}]);
                    [Gpath, Gfile] = xASL_fileparts(dest_file);
                    dest_BIDS = fullfile(Gpath, [Gfile BIDSext{iB}]);
                    if exist(temp_BIDS, 'file')
                        xASL_Move(temp_BIDS, dest_BIDS, parms.Overwrite, parms.Verbose);
                    end
                end
            end
        end
    catch err
        rmdir(temp_dir,'s');
        rethrow(err)
    end

    %% cleanup temp
    rmdir(temp_dir,'s');

    %% optionally return the used input file
    if nargout>1
        usedinput = inpath; % could be the original input or the dicom file used
    end
end
