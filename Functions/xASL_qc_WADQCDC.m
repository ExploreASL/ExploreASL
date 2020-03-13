function xASL_qc_WADQCDC(x, iSubject, ScanType)
%xASL_qc_WADQCDC Runs WAD-QC specific Python script to incorporate QC parameters in DICOM
%
% FORMAT: xASL_qc_WADQCDC(x, iSubject[, ScanType])
%
% INPUT:
%   x           - structure containing fields with all information required to run this submodule (REQUIRED)
%   iSubject    - index of current subject (REQUIRED)
%   ScanType    - ScanType to process (OPTIONAL, DEFAULT=loop over all ScanTypes: {'Structural', 'ASL', 'func', 'dwi'})
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This QC function runs WAD-QC specific Python script to zip QC information &
%              incorporate this into a DICOM field for analysis on the
%              WAD-QC server, by the following:
%              a) Define QCDC script: this is the Python script written by
%                 Gaspare, edited by Joost Kuijer & copied to the EPAD
%                 CustomScripts folder of ExploreASL
%              b) Python installation location is checked, with several known
%                 locations, for several servers. If we cannot find it,
%                 the QCDC is not ran
%              c) Previous QCDC results are cleaned. QCDC stores all its
%                 results in a separate folder (Something like 2 layers up from the current
%                 folder, here referred to as QCDCDir = [x.D.ROOT 'qcdc_output'])
%                 from these result files, only the filled DICOM file is
%                 interesting, all the rest are copies of the QC results
%                 that we embedded into the DICOM
%              d) Run QCDC (if Python installation detected)
%                 The following files need to be set as executable:
%                     ('QCDC', 'src', 'qc_data_collector.py')
%                     ('QCDC', 'src', 'bash', 'create_dcm_only_wadqc.sh')
%                     ('QCDC', 'src', 'bash', 'sendwadqc.sh')
%              e) Clean up new QCDC results (as above) & move the filled
%                 DICOM to ['qcdc_' DummyFile] within the current ScanType
%                 folder
%              f) Sending the DICOM to the WAD-QC server using storescu
%
% EXAMPLE: xASL_qc_WADQCDC(x, 10, 'ASL');
% __________________________________
% Copyright (C) 2015-2019 ExploreASL



%% Loop over Scantypes
ScanTypes = {'Structural', 'ASL', 'func', 'dwi'};
SubPath = {x.SUBJECTS{iSubject}, fullfile(x.SUBJECTS{iSubject}, 'ASL_1'), fullfile(x.SUBJECTS{iSubject}, 'func'), fullfile(x.SUBJECTS{iSubject}, 'dwi')};
NIfTIname = {'T1', 'ASL4D', 'func', 'dwi'};

if nargin>2 && ~isempty(ScanType)
    % if we have specified Current ScanType (i.e. when running within a
    % specific module)
    CurrInd = find(strcmp(ScanTypes,ScanType));
    ScanTypes = {ScanTypes{CurrInd}};
    SubPath = {SubPath{CurrInd}};
    NIfTIname = {NIfTIname{CurrInd}};
end % otherwise, simply iterate over all scantypes below

for iScanType=1:length(ScanTypes)

    WAD_Path = fullfile(x.D.ROOT, SubPath{iScanType}, ['qcdc_' x.SUBJECTS{iSubject} '_' ScanTypes{iScanType} '.json']);
    DummyDir = fullfile(x.D.ROOT, SubPath{iScanType});
    DummyFile = xASL_adm_GetFileList(DummyDir, ['^DummyDicom_' NIfTIname{iScanType} '.*.dcm$'],'List',[0 Inf]);

    if ~exist(WAD_Path,'file')
        warning([WAD_Path ' missing, skipping qcdc']);
        continue;
    elseif isempty(DummyFile)
        warning('DummyDicom missing, skipping qcdc');
        continue;
    elseif length(DummyFile)>1
        warning('Multiple DummyDicoms found, using first for qcdc');
    end
    DummyFile = DummyFile{1};
    DummyPath = fullfile(DummyDir, DummyFile);

%         %% Check if the AcquisitionDate parameter exists
%         Info = xASL_io_ReadTheDicom(false, DummyPath); % try to use DCM_TK
%
%         if isfield(Info,'AcquisitionDate') && ~isempty(Info.AcquisitionDate)
%             ExistDCMAcqDateField = true;
%         else
%             ExistDCMAcqDateField = false;
%         end
%         if isfield(Info,'StudyDate') && ~isempty(Info.StudyDate)
%             ExistDCMStudyDateField = true;
%         else
%             ExistDCMStudyDateField = false;
%         end

    %% a) Define Gaspare's script (QCDC)
    QCDC_Path = fullfile(x.MyPath, 'CustomScripts', 'EPAD', 'QCDC', 'src', 'qc_data_collector.py');
    QCDC_sh1 = fullfile(x.MyPath, 'CustomScripts', 'EPAD', 'QCDC', 'src', 'bash', 'create_dcm_only_wadqc.sh');
    QCDC_sh2 = fullfile(x.MyPath, 'CustomScripts', 'EPAD', 'QCDC', 'src', 'bash', 'sendwadqc.sh');

    %% b) Check for Python installation
    [StatusPython, dummy] = system('ls -d /data/usr/local/anaconda2/bin'); % Victory VUmc cloud
    [StatusPython2, dummy] = system('wsl ls -d /usr/src/Python-2.7.16'); % Victory VUmc cloud
    [StatusPython3, dummy] = system('ls -d /usr/bin/python'); % AMC flux

    if StatusPython==0
        PythonPath = '/data/usr/local/anaconda2/bin/python';
        RunPython = true;
    elseif StatusPython2==0
        PythonPath = 'wsl /usr/src/Python-2.7.16/python';
        RunPython = true;
    elseif StatusPython3==0
        PythonPath = '/usr/bin/python2.7';
        RunPython = true;
    else
        RunPython = false;
    end

    %% c) CleanUp previous QCDC results
    QCDCDir = [x.D.ROOT 'qcdc_output_qcdc_' x.SUBJECTS{iSubject} '_' ScanTypes{iScanType}];
    xASL_adm_DeleteFileList(QCDCDir, '.*', true, [0 Inf]);

    %% d) Run QCDC
    if RunPython
        fprintf('Incorporating QC results into dummy dicom...\n');
        cd(x.D.ROOT);
        system(['chmod +x ' QCDC_Path]);
        system(['chmod +x ' QCDC_sh1]);
        system(['chmod +x ' QCDC_sh2]);
        result = system([PythonPath ' ' xASL_adm_UnixPath(QCDC_Path) ' ' xASL_adm_UnixPath(WAD_Path) ' ' xASL_adm_UnixPath(x.D.ROOT)]);

        if result~=0
            % throw error
            fprintf('Couldnt incorporate QC results in dicom, QCDC failed\n');
        end

        %% e) Clean up QCDC output results
        TempPath = fullfile(QCDCDir, DummyFile);
        NewPath = fullfile(x.D.ROOT, SubPath{iScanType}, ['qcdc_' DummyFile]);
        if exist(TempPath, 'file')
            xASL_Move(TempPath, NewPath, true);
            fprintf('QCDC has succesfully incorporated QC into dummy dicom\n');
            IsSuccess = true;
        else
            fprintf('Imported WAD-QC DICOM missing, QCDC failed\n');
            IsSuccess = false;
        end

        xASL_delete(QCDCDir, true);

        %% f) Try sending the DICOM to WAD-QC server
        if IsSuccess
            cd(x.D.ROOT);
            [results, results2] = system(['storescu wad-qc 11112 ' NewPath]);
            if results==0 && isempty(results2)
                fprintf('DICOM sent to WADQC server\n');
            else
                fprintf('Couldnt send DICOM to WADQC server:\n');
                fprintf(results2);
            end
        end
    end
end


end
