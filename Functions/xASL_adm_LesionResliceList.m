function [INname, OUTname] = xASL_wrp_LesionResliceList(x,bLesion_T1,bLesion_FLAIR,bROI_T1,bROI_FLAIR)
%xASL_wrp_LesionResliceList Creates list of structural image paths to
% reslice.
%
% FORMAT:       [INname, OUTname] = xASL_wrp_LesionResliceList(x,bLesion_T1,bLesion_FLAIR,bROI_T1,bROI_FLAIR)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Creates list of structural image paths to reslice.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL

% Do everything unless specified otherwise
if nargin<2 || isempty(bLesion_T1)
	bLesion_T1 = 1;
end

if nargin<3 || isempty(bLesion_FLAIR)
	bLesion_FLAIR = 1;
end

if nargin<4 || isempty(bROI_T1)
	bROI_T1 = 1;
end

if nargin<5 || isempty(bROI_FLAIR)
	bROI_FLAIR = 1;
end

INname      = '';
OUTname     = '';

% Load those lists that are asked
if bLesion_T1
	Lesion_T1_list      = xASL_adm_GetFileList(x.dir.SUBJECTDIR, ['^Lesion_' x.P.STRUCT '_\d*\.(nii|nii\.gz)$'], 'FPList', [0 Inf]);
else
	Lesion_T1_list = '';
end

if bLesion_FLAIR
	Lesion_FLAIR_list   = xASL_adm_GetFileList(x.dir.SUBJECTDIR, ['^Lesion_' x.P.FLAIR '_\d*\.(nii|nii\.gz)$'], 'FPList', [0 Inf]);
else
	Lesion_FLAIR_list = '';
end

if bROI_T1
	ROI_T1_list         = xASL_adm_GetFileList(x.dir.SUBJECTDIR, ['^ROI_' x.P.STRUCT '_\d*\.(nii|nii\.gz)$'], 'FPList', [0 Inf]);
else
	ROI_T1_list = '';
end

if bROI_FLAIR
	ROI_FLAIR_list      = xASL_adm_GetFileList(x.dir.SUBJECTDIR, ['^ROI_' x.P.FLAIR '_\d*\.(nii|nii\.gz)$'], 'FPList', [0 Inf]);
else
	ROI_FLAIR_list = '';
end

if  exist('Lesion_T1_list', 'var')
    for iS=1:length(Lesion_T1_list)
        INname{end+1}       = Lesion_T1_list{iS};

        [~, Ffile, ~]       = xASL_fileparts(Lesion_T1_list{iS});
        ModLesion           = x.P.STRUCT;
        nLesion             = num2str(Ffile(9+length(ModLesion):end));

        OUTname{end+1}      = fullfile(x.D.PopDir, ['rLesion_' ModLesion '_' nLesion '_' x.P.SubjectID '.nii']);
    end
end

if  exist('Lesion_FLAIR_list', 'var')    
    for iS=1:length(Lesion_FLAIR_list)
        INname{end+1}       = Lesion_FLAIR_list{iS};

        [~, Ffile, ~]       = xASL_fileparts(Lesion_FLAIR_list{iS});
        ModLesion           = x.P.FLAIR;
        nLesion             = num2str(Ffile(9+length(ModLesion):end));

        OUTname{end+1}      = fullfile(x.D.PopDir, ['rLesion_' ModLesion '_' nLesion '_' x.P.SubjectID '.nii']);
    end    
end

if  exist('ROI_T1_list', 'var')
    for iS=1:length(ROI_T1_list)
        INname{end+1}       = ROI_T1_list{iS};

        [~, Ffile, ~]       = xASL_fileparts(ROI_T1_list{iS});
        ModROI              = x.P.STRUCT;
        nROI                = num2str(Ffile(6+length(ModROI):end));

        OUTname{end+1}      = fullfile(x.D.PopDir, ['rROI_' ModROI '_' nROI '_' x.P.SubjectID '.nii']);
    end
end

if  exist('ROI_FLAIR_list', 'var')    
    for iS=1:length(ROI_FLAIR_list)
        INname{end+1}       = ROI_FLAIR_list{iS};

        [~, Ffile, ~]       = xASL_fileparts(ROI_FLAIR_list{iS});
        ModROI              = x.P.FLAIR;
        nROI                = num2str(Ffile(6+length(ModROI):end));

        OUTname{end+1}      = fullfile(x.D.PopDir, ['rROI_' ModROI '_' nROI '_' x.P.SubjectID '.nii']);
    end    
end


end