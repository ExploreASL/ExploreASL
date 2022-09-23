function [niifiles, ScanNameOut, usedinput, msg] = xASL_io_dcm2nii(inpath, destdir, series_name, imPar, myPath)
%xASL_io_dcm2nii Convert DICOM NIfTI/BIDS format using the dcm2nii command line utility.
% (http://www.nitrc.org/projects/mricron)
%
% FORMAT:       [niifiles, ScanNameOut, usedinput, msg] = xASL_io_dcm2nii(inpath, destdir, series_name, imPar, myPath)
% 
% INPUT:   
%      inpath          path to dicom folder, dicom file, PAR-file or REC-file. In case of a dicom folder, 'DcmExt' will be used to
%                      filter the folder contents. In case of a dicom file, all dicom files in the same folder will be used. (REQUIRED)
%      destdir         target destination directory (REQUIRED)
%	   series_name     target name of the NIfTI file (REQUIRED)
%      imPar           imPar from x.modules.import.imPar (OPTIONAL, DEFAULT = [])
%                      'imPar.dcm2nii_version' (DEFAULT 20220720)
%                      'imPar.bVerbose'        set to true or false to switch terminal feedback (DEFAULT false)
%                      'imPar.bOverwrite'      set to true or false to overwrite existing files (DEFAULT false)
%                      'imPar.dcmExtFilter' regular expression used to find dicom files (DEFAULT '^.+\.dcm$')
%      myPath          myPath to the xASL folder from x.opts.MyPath (OPTIONAL, DEFAULT = pwd)
%
%               The first DICOM file will be used if inpath is a directory.
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Convert DICOM NIfTI/BIDS format using the dcm2nii command line utility.
%
% 1. Initial settings
% 2. Parse parameters
% 3. Locate dcm2nii executable
% 4. Check if we are reading a DICOM folder
% 5. Set dcm2niiX initialization loading
% 6. Check for existing targets
% 7. Create temporary subfolder for converting
% 8. Run dcm2nii and move files to final destination using specified series name
% 9. Cleanup temp
% 10. Optionally return the used input file
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2022 ExploreASL

    %% 1. Initial settings

	if nargin < 3
		error('Require at least three input parameters inpath, destdir, series_name');
	end
	
	if nargin < 4 || isempty(imPar)
		imPar = struct;
	end
	
	if ~isfield(imPar,'dcm2nii_version') || isempty(imPar.dcm2nii_version)
		imPar.dcm2nii_version = '20220720';
	end
	
	if ~isfield(imPar, 'dcmExtFilter') || isempty(imPar.dcmExtFilter)
		imPar.dcmExtFilter = '^.+\.dcm$';
	end
	
	if ~isfield(imPar, 'bVerbose') || isempty(imPar.bVerbose)
		imPar.bVerbose = false;
	end
	
	if ~isfield(imPar, 'bOverwrite') || isempty(imPar.bOverwrite)
		imPar.bOverwrite = false;
	end
	
	if nargin < 5 || isempty(myPath)
		myPath = pwd;
	end
	
    niifiles = {};
    msg = [];

    %% 3. Locate dcm2nii executable
	
	if ismac && str2num(imPar.dcm2nii_version(1:4))<2014
        imPar.dcm2nii_version = '20220720'; % mac is incompatible with older versions
	end

    mricron_path = fullfile(myPath,'External','MRIcron');
	mricron_version_path = fullfile(mricron_path, imPar.dcm2nii_version);
	switch computer()
		case {'PCWIN', 'PCWIN64'}
			ExePath = fullfile(mricron_version_path,'dcm2niix.exe');
		case 'GLNX86'
			ExePath = fullfile(mricron_version_path,'dcm2niix-lx32');
		case 'GLNXA64'
			ExePath = fullfile(mricron_version_path,'dcm2niix-lx64');
		case {'MACI','MACI64'}
			ExePath = fullfile(mricron_version_path,'dcm2niix-osx');
		otherwise
			error('Unknown computer architecture: %s',computer());
	end
    [stat, attr] = fileattrib(ExePath);
    if ~stat
        error('dcm2niix application not found: %s',ExePath);
    end
    if ~attr.UserExecute
        % try fixing
        try
            system(['chmod 777 ' ExePath]);
        catch
            disp(attr)
            error('dcm2niix application is not executable!');
        end
    end
	
    %% 4. Check if we are reading a DICOM folder
    if exist(inpath,'dir')
        % Check if there are any DICOM files
        dicom_files = xASL_adm_GetFileList(inpath, imPar.dcmExtFilter,'List');
        if isempty(dicom_files)
            % This error will be catched in wrapper function, verbose & error catching will be dealt with there
            error(['No dicom files match ' imPar.dcmExtFilter ' in ' inpath ', skipping...']);
        end

        %% sort the dicom filenames because it's nice to use the first (but not required!)
        dicom_files = sort(dicom_files);
        inpath = fullfile(inpath, dicom_files{1});
        bScanDicomFolder = true;
    else
        bScanDicomFolder = false;
    end

	if bScanDicomFolder
		IniPath = fullfile(mricron_path,'dcm2nii-dicom.ini');
	else
		IniPath = fullfile(mricron_path,'dcm2nii-parrec.ini');
	end
	
	%% 5. Set dcm2niiX initialization loading
    if str2num(imPar.dcm2nii_version(1:4))<2014
        dcm2nii_args = sprintf('-b "%s"', IniPath);
    else
        dcm2nii_args = '-f "%f_%p_%t_%s_%r"';
    end

    %% 6. Check for existing targets
    bSingleExists = xASL_exist(fullfile(destdir, [series_name '.nii']),'file');
    bMultiExists = xASL_exist(fullfile(destdir, [series_name '_1.nii']),'file');
    if ~imPar.bOverwrite && (bSingleExists || bMultiExists)
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

        if imPar.bVerbose; fprintf('%s\n', [ScanNameOut ' already existed, skipping...']); end
        return;
    end


    %% 7. Create temporary subfolder for converting
    temp_dir = fullfile(destdir, ['dcm2nii_temp_' series_name]);
    
    % Ignore previously stored temp files that are still here because of a code crash e.g.
    if exist(temp_dir,'dir')
        xASL_adm_DeleteFileList(temp_dir,'.*', true, [0 Inf]);
    end
    xASL_adm_CreateDir(temp_dir);

    %% 8. Run dcm2nii and move files to final destination using specified series name
    
    % Use catch to remove temporary folder on error
    try
        % OS check
        if ispc()
            quote = '"';
        else
            quote = '''';
        end

        % Command line string
        if ismac()
            cmdline = [ExePath ' -o ' temp_dir ' ' inpath];
        else
            cmdline = [quote ExePath quote ' ' dcm2nii_args ' -o ' quote temp_dir quote ' ' quote inpath quote];
        end

        % User feedback if verbose
        if imPar.bVerbose
            fprintf('executing: [%s]\n',cmdline);
            [status, msg] = system(cmdline, '-echo');
            separatorline = xASL_adm_BreakString('',[],[],1,0);
            fprintf('%s\n%s%s\nstatus: %d\n',separatorline,msg,separatorline,status);
        else
            [~, msg] = system(cmdline);
        end

        %% Move/Rename NIfTIs to final destination
        niiEntries = xASL_adm_GetFileList(temp_dir, '.*\.nii$', 'FPListRec', [0 Inf]);
		
		% Check if NIfTIs do not exist
		if isempty(niiEntries)
            
			% Try to read files with LS command
			pathEntriesAlt = ls(temp_dir);
			
            % Not yet properly tested on Windows, as the error hasn't occured there with the same data as on Linux
            if ~isempty(pathEntriesAlt)
				if ispc
					warning('No NIfTIs found! Either illegal characters in filenames or other code defect...');
				end
                % Then break to separate files according to character 32
                indNewLine = find(pathEntriesAlt == 10);

                %indNewLine = [indNewLine(:);length(pathEntriesAlt)+1];
                indexStart = 1;
                niiEntriesAlt = '';
                for iEntry = 1:length(indNewLine)
                    % Skip the new line cr (10 32) at the end
                    niiEntriesAlt{iEntry} = pathEntriesAlt(indexStart:(indNewLine(iEntry)-1));
                    indexStart = indNewLine(iEntry)+1;
                end

                for iEntry = 1:length(niiEntriesAlt)
                    % Remove leading and trailing apostrophes and 10s for
                    if niiEntriesAlt{iEntry}(1) == 32
                        niiEntriesAlt{iEntry} = niiEntriesAlt{iEntry}(2:end);
                    end
                    if niiEntriesAlt{iEntry}(end) == 32
                        niiEntriesAlt{iEntry} = niiEntriesAlt{iEntry}(1:(end-1));
                    end

                    if niiEntriesAlt{iEntry}(1) == 39
                        niiEntriesAlt{iEntry} = niiEntriesAlt{iEntry}(2:end);
                    end
                    if niiEntriesAlt{iEntry}(end) == 39
                        niiEntriesAlt{iEntry} = niiEntriesAlt{iEntry}(1:(end-1));
                    end

                    % Look for further '$' and '' as starting and ending patterns
    				indexStart = strfind(niiEntriesAlt{iEnt}, string('''$'''));
    				indexEnd = strfind(niiEntriesAlt{iEnt}, string(''''''));

                    % If both patterns exist, then remove its contents (insides)
                    
                    if ~isempty(indexStart) && ~isempty(indexEnd)
                        % Within that count the number of special characters by backslashes
                        numSlashes = sum(niiEntriesAlt{iEntry}((indexStart+3):(indexEnd-1)) == '\');

                        % Replace by ? and move files to a normal name
                        fTempNiiSubstitute = [niiEntriesAlt{iEntry}(1:(indexStart-1)), repmat('?',[1 numSlashes]), niiEntriesAlt{iEntry}((indexEnd+2):end)];
                        fTempNiiFixed = niiEntriesAlt{iEntry}([1:(indexStart-1),(indexEnd+2):end]);
                        xASL_SysMove(fTempNiiSubstitute,fTempNiiFixed,[],false);
                    end
                end
                % Read the files again
                niiEntries = xASL_adm_GetFileList(temp_dir, '.*\.nii$', 'FPListRec', [0 Inf]);
            end
		end
		
        if isempty(niiEntries)
            error('Empty output dcm2niix');
        else
            [niiPaths, niiNames, ~] = cellfun(@(y) xASL_fileparts(y), niiEntries, 'UniformOutput',false);
        end

        % Number of NIfTIs
        nNifties = length(niiNames);
        if nNifties==0
            error('Conversion failed: No nifti output files available');
		else
			vectorKeep = 1:nNifties;
            
            % Fix NIFTI file names
            [niiEntries, niiNames] = xASL_io_dcm2nii_FixFormat(niiEntries, niiPaths, niiNames);

            % Check if dcm2niix added some contrast information at the end of the FileName
            if length(vectorKeep)==1
                ContrastAdd{1} = [];
            else
                % first compare the filenames (for their equivalent length)
                EquiL = min(cellfun(@(y) length(y), niiNames));

                NN=1;
                for iN1=vectorKeep
                    for iN2=vectorKeep
                        CmpString(NN,:) = double(niiNames{iN1}(1:EquiL)~=niiNames{iN2}(1:EquiL));
                        NN=NN+1;
                    end
                end
                IndCommon = find(sum(CmpString,1)==size(CmpString,1)); % common filename differences (i.e. within equivalent filename length)
                for iN=vectorKeep
                    ContrastAdd{iN} = niiNames{iN}(IndCommon);
                    if length(niiNames{iN})>EquiL
                        ContrastAdd{iN} = [ContrastAdd{iN} niiNames{iN}(EquiL+1:end)];
                    end
                end
            end

            %% Iterate over volumes
            for iVolume=sort(vectorKeep)
                fTempNii = niiEntries{iVolume};
				
                % Make BIDS compatible dest_file here
                DestFileName = series_name;

                %[Ind1, Ind2] = regexp(niiNames{iVolume},'_run-\d*_','start','end');
                %if ~isempty(Ind1) && ~isempty(Ind2) % add run_index suffix (if it is there)
                %    RunName = niiNames{iVolume}(Ind1+1:Ind2-1);
                %    DestFileName = [DestFileName '_' RunName];
                %end
                %if ~isempty(ContrastAdd{iVolume}) % add contrast suffix (if there is something that indicates it)
                %    DestFileName = [DestFileName '_' ContrastAdd{iVolume}];
                %end

                if iVolume==min(vectorKeep)
                    ScanNameOut = DestFileName;
                end
                
                if length(vectorKeep)>1 % add iVolume suffix for SeriesNumber (if there are multiple)
                    % Obtain the SeriesNumber from JSON
                    [Gpath, Gfile] = xASL_fileparts(fTempNii);
                    tempJSON = fullfile(Gpath,[Gfile '.json']);
                    if exist(tempJSON,'file')
                        tempJSON = spm_jsonread(tempJSON);
                        if isfield(tempJSON,'SeriesNumber')
                            DestFileName = [DestFileName '_' num2str(tempJSON.SeriesNumber)];
                        end
                    else
                        warning('JSON sidecar missing after dcm2niix conversion');
                    end
                end
                
                % Fallback
                niiInstanceNumber = [];
                
                % Check for files that are formatted like ..._InstanceNumber_eNumber instead of ..._InstanceNumber
                if isempty(niiInstanceNumber)
                    expression = '_(\d+)_e.+$';
                    [~, fileName, ~] = xASL_fileparts(fTempNii);
                    startIndex = regexp(fileName,expression);
                    if ~isempty(startIndex)
                        niiInstanceNumber = fileName(startIndex+1:end);
                        niiInstanceNumber = niiInstanceNumber(1:strfind(niiInstanceNumber,'_')-1);
                    end
                end
				
                % Determine InstanceNumbers from fileNames
                if isempty(niiInstanceNumber)
                    expression = '_(\d+)$'; % Get last number after last _ symbol
                    [~, fileName, ~] = xASL_fileparts(fTempNii);
                    startIndex = regexp(fileName,expression);
                    niiInstanceNumber = fileName(startIndex+1:end);
                end

                if length(vectorKeep)>1 % add iVolume suffix (if there are multiple)
                    try
                        DestFileName = [DestFileName '_' niiInstanceNumber];
                    catch
                        fprintf('Something went wrong while trying to get the InstanceNumber...\n');
                        DestFileName = [DestFileName '_' int2str(iVolume)];
                    end
                end
                
                % DestFileName = xASL_adm_CorrectName(DestFileName); %
                % TAKING THIS OUT HERE, IF WE WANT THIS THIS NEEDS TO BE
                % EQUAL FOR THE dcm2niiX OUTPUT AND THE JSON THAT WE CREATE
                % WITH DCMTK
                % ALSO BIDS REQUIRES TO KEEP '-' IN
                
                dest_file = fullfile(destdir,[DestFileName '.nii']);
				
				% Check for suspicious illegal characters
                indIchar = find((fTempNii < 32) | (fTempNii > 126));
                if ~isempty(indIchar)
                    % If these characters are present in the file name, we
                    % will not be able to move it in the following step
                   
                    % We therefore replace the illegal characters by ? to
                    % allow moving
                    fTempNiiSubstitute = fTempNii;
                    fTempNiiSubstitute(indIchar) = '?';
                    fTempNiiFixed = fTempNii;
                    fTempNiiFixed(indIchar) = '_';
                    
                    % And we then move to a file with a similar name, but
                    % illegal characters replaced by a
                    xASL_SysMove(fTempNiiSubstitute,fTempNiiFixed,[],false);
                    
                    % From then on, work with the corrected file name
                    fTempNii = fTempNiiFixed;
                    
                    % Apply the same to other BIDS files
                    BIDSext = {'.json' '.bval' '.bvec'};
                    for iB=1:length(BIDSext)
                        [Gpath, Gfile] = xASL_fileparts(fTempNiiSubstitute);
                        temp_BIDS = fullfile(Gpath, [Gfile BIDSext{iB}]);
                        [Gpath, Gfile] = xASL_fileparts(fTempNiiFixed);
                        dest_BIDS = fullfile(Gpath, [Gfile BIDSext{iB}]);
                        xASL_SysMove(temp_BIDS, dest_BIDS, [],false);
                    end
                end
				
                xASL_Move(fTempNii, dest_file, imPar.bOverwrite, imPar.bVerbose);
                niifiles{end+1} = dest_file; %#ok<AGROW>

                % Do the same for the BIDS files
                BIDSext = {'.json' '.bval' '.bvec'};
                for iB=1:length(BIDSext)
                    [Gpath, Gfile] = xASL_fileparts(fTempNii);
                    temp_BIDS = fullfile(Gpath, [Gfile BIDSext{iB}]);
                    [Gpath, Gfile] = xASL_fileparts(dest_file);
                    dest_BIDS = fullfile(Gpath, [Gfile BIDSext{iB}]);
                    % Move JSON file
                    if xASL_exist(temp_BIDS, 'file')
                        xASL_Move(temp_BIDS, dest_BIDS, imPar.bOverwrite, imPar.bVerbose);
                        % Add InstanceNumber if possible
                        if ~isempty(niiInstanceNumber)
                            tmpJSON = spm_jsonread(dest_BIDS);
                            tmpJSON.InstanceNumber = niiInstanceNumber;
                            spm_jsonwrite(dest_BIDS,tmpJSON);
                        end
                    end
                end
            end
        end
    catch err
        rmdir(temp_dir,'s');
        rethrow(err)
    end

    %% 9. Cleanup temp
    rmdir(temp_dir,'s');

    %% 10. Optionally return the used input file
    if nargout>1
        usedinput = inpath; % could be the original input or the dicom file used
    end
end


%% Fix NIfTI and JSON file names (uneven length etc.)
function [niiEntriesDest, niiNamesDest] = xASL_io_dcm2nii_FixFormat(niiEntries, niiPaths, niiNames)

    % Fallback
    instanceNumber = '';
    echoNumber = '';
    niiNamesDest = niiNames;
    niiEntriesDest = niiEntries;

    % Possible expressions
    expressionA = '_(\d+)_e.$';
    expressionB = '_(\d+)$';

    % Iterate over files
    for iFile = 1:size(niiNames,1)
        
        % Check possible formats (try B if A does not work)
        startIndexA = regexp(niiNames{iFile},expressionA);
        startIndexB = regexp(niiNames{iFile},expressionB);
        
        % Try different formats
        if ~isempty(startIndexA)
            elements = strsplit(niiNames{iFile}(startIndexA+1:end),'_');
            instanceNumber = elements{1};
            echoNumber = elements{2}(2:end);
        elseif ~isempty(startIndexB)
            instanceNumber = niiNames{iFile}(startIndexB+1:end);
        end
        
        % Use zero padding
        instanceNumber = sprintf('%05s', instanceNumber);
        echoNumber = sprintf('%05s', echoNumber);
        
        % Rename file
        if ~isempty(startIndexA)
            niiNamesDest{iFile} = [niiNames{iFile}(1:startIndexA) instanceNumber '_e' echoNumber];
        elseif ~isempty(startIndexB)
            niiNamesDest{iFile} = [niiNames{iFile}(1:startIndexB) instanceNumber '_e' echoNumber];
        else
            niiNamesDest{iFile} = [niiNames{iFile} instanceNumber '_e' echoNumber];
        end
        
        % Actual rename
        niiEntriesDest{iFile} = fullfile(niiPaths{iFile},[niiNamesDest{iFile} '.nii']);
        xASL_Move(fullfile(niiPaths{iFile},[niiNames{iFile} '.nii']),fullfile(niiPaths{iFile},[niiNamesDest{iFile} '.nii']));
        xASL_Move(fullfile(niiPaths{iFile},[niiNames{iFile} '.json']),fullfile(niiPaths{iFile},[niiNamesDest{iFile} '.json']));
        
    end

end
