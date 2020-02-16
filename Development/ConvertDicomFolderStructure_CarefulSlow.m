function ConvertDicomFolderStructure_CarefulSlow(ROOT, bUseDCMTK, bVerbose)
%ConvertDicomFolderStructure_CarefulSlow Script to put dicom names in directories according to their ProtocolName/SeriesDescrption

if nargin<2 || isempty(bUseDCMTK)
    bUseDCMTK = true; % faster
end
if nargin<3 || isempty(bVerbose)
    bVerbose = true;
end

    Flist   = xASL_adm_GetFileList(ROOT,'^.*.(?!(xlsx|ini))$','FPListRec',[0 Inf]);

    for iL=1:length(Flist)
        if bVerbose; xASL_TrackProgress(iL,length(Flist)); end
        clear tDcm Fname NewDir NewFile Ppath Pfile Pext
        
        try
            tDcm = ReadTheDicom(bUseDCMTK, Flist{iL});

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

            [Ppath Pname Pext]  = fileparts(Flist{iL});
            NewFile             = fullfile(NewDir,[Pname '.dcm']);

            if  exist(Flist{iL},'file') && ~exist(NewFile,'file')
                xASL_Move(Flist{iL},NewFile);
            end
        end
    end

    Dlist   = xASL_adm_GetFsList(ROOT,'^.*$',1,0,0,[0 Inf]);
    for iD=3:length(Dlist)
        if  isempty(xASL_adm_GetFileList(fullfile(ROOT,Dlist{iD}),'^.*','FPListRec',[0 Inf]))
            try
                rmdir(fullfile(ROOT,Dlist{iD}));
            end
        end
    end
end



%% =============================================================================================================================                
%% =============================================================================================================================
function [Info] = ReadTheDicom(bUseDCMTK, DicomPath)
%ReadTheDicom Wrapper around DICOMTK & dicominfo for reading the dicom fields

warning('off','images:dicominfo:fileVRDoesNotMatchDictionary');

try
    if bUseDCMTK
        Info = xASL_io_DcmtkRead(DicomPath);

        if isempty(Info.EchoTime) || isempty(Info.RepetitionTime) || isempty(Info.ImageType)
            Info = dicominfo(DicomPath);
            if ~isfield(Info,'EchoTime') || ~isfield(Info,'RepetitionTime') || ~isfield(Info,'ImageType')
                error('Both DCMTK output & dicominfo were incomplete');
            else
                error('DCMTK output was incomplete, but dicominfo was complete');
            end
        end
    else
        Info = dicominfo(DicomPath);
    end
catch ME
    warning(ME.message);
end

end