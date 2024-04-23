function ConvertDicomFolderStructure_CarefulSlow(ROOT, bUseDCMTK, bVerbose)
%ConvertDicomFolderStructure_CarefulSlow Script to put dicom names in directories according to their ProtocolName/SeriesDescrption
% Also add .dcm extensions

if nargin<2 || isempty(bUseDCMTK)
    bUseDCMTK = true; % faster
end
if nargin<3 || isempty(bVerbose)
    bVerbose = true;
end

    Flist   = xASL_adm_GetFileList(ROOT,'^.*.(?!(xlsx|ini))$','FPListRec',[0 Inf]);

    for iL=1:length(Flist)
        if bVerbose; xASL_TrackProgress(iL, length(Flist)); end
        clear tDcm Fname NewDir NewFile Ppath Pfile Pext
        
        try
			tDcm = xASL_io_DcmtkRead(Flist{iL}, false, bUseDCMTK);

			if isempty(tDcm.EchoTime) || isempty(tDcm.RepetitionTime) || isempty(ItDcm.ImageType)
				warning(['Incomplete DICOM header: ' DicomPath]);
			end




            try
                if ~strcmp(tDcm.ProtocolName,tDcm.SeriesDescription)
                    Fname = [tDcm.ProtocolName '_' tDcm.SeriesDescription];
                else
                    Fname = tDcm.ProtocolName;
                end
            catch
                Fname = tDcm.ProtocolName;
            end
            
            Fname = xASL_adm_CorrectName(Fname);
            NewDir = fullfile(ROOT,Fname);
            xASL_adm_CreateDir(NewDir);

            [~, Pname, Pext] = fileparts(Flist{iL});
            if strcmpi(Pext, '.ima') || strcmpi(Pext, '.dcm')
                NewFile = Flist{iL}; % don't change the filename
            else
                % in the case of a period (.) in the filename, the extension may be
                % incorrectly detected and we need to append .dcm
                NewFile = fullfile(NewDir,[Pname Pext '.dcm']);
            end

            if ~strcmp(Flist{iL},NewFile) && exist(Flist{iL},'file') && ~exist(NewFile,'file')
                xASL_Move(Flist{iL}, NewFile);
            end
        end
    end

    Dlist = xASL_adm_GetFsList(ROOT,'^.*$',1,0,0,[0 Inf]);
    for iD=3:length(Dlist)
        if  isempty(xASL_adm_GetFileList(fullfile(ROOT,Dlist{iD}),'^.*','FPListRec',[0 Inf]))
            try
                rmdir(fullfile(ROOT,Dlist{iD}));
            end
        end
    end
end
