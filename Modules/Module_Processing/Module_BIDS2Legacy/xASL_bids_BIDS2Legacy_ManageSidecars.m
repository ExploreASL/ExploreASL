function [bidsPar, pathOrig, pathDest, TypeIs] = xASL_bids_BIDS2Legacy_ManageSidecars(bidsPar, pathOrig, pathDest, TypeIs)
%xASL_bids_BIDS2Legacy_ManageSidecars Manage JSON sidecars for BIDS2Legacy conversion.
%
% FORMAT:     [bidsPar, pathOrig, pathDest, TypeIs] = xASL_bids_BIDS2Legacy_ManageSidecars(bidsPar, pathOrig, pathDest, TypeIs)
% 
% INPUT:      bidsPar  - BIDS par struct (STRUCT, REQUIRED)
%             pathOrig - Origin path (STRING, REQUIRED)
%             pathDest - Destination path (STRING, REQUIRED)
%             TypeIs   - Type (REQUIRED)
%   
% OUTPUT:     bidsPar  - BIDS par struct (STRUCT, REQUIRED)
%             pathOrig - Origin path (STRING, REQUIRED)
%             pathDest - Destination path (STRING, REQUIRED)
%             TypeIs   - Type (REQUIRED)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Manage JSON sidecars for BIDS2Legacy conversion.
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     [bidsPar, pathOrig, pathDest, TypeIs] = xASL_bids_BIDS2Legacy_ManageSidecars(bidsPar, pathOrig, pathDest, TypeIs);
% __________________________________
% Copyright 2015-2021 ExploreASL


    iCount = 1;
    for iCar=1:length(bidsPar.sidecarName)
        [Fpath, Ffile] = xASL_fileparts(pathOrig{1});

        if ~bidsPar.sidecarSuffixType(iCar)
            Ffile = Ffile(1:end-length(TypeIs)-1);
        end
        TempSidecar = fullfile(Fpath, [Ffile bidsPar.sidecarName{iCar}]);

        if ~strcmp(bidsPar.sidecarTypeSpecific{iCar}, 'no') && ~strcmp(bidsPar.sidecarTypeSpecific{iCar}, TypeIs)
            % skip this sidecar (e.g. some asl-specific sidecars for non-asl NIfTIs)
        elseif ~exist(TempSidecar, 'file') && bidsPar.sidecarRequired(iCar)
            warning([TempSidecar ' missing']);
        elseif exist(TempSidecar, 'file')
            pathOrig{iCount+1} = TempSidecar;

            [Fpath, Ffile] = xASL_fileparts(pathDest{1});
			switch (bidsPar.sidecarName{iCar})
				case '.json' % Copy side-car name unchanged
					pathDest{iCount+1} = fullfile(Fpath, [Ffile bidsPar.sidecarName{iCar}]);
				case '_aslcontext.tsv'
					pathDest{iCount+1} = fullfile(Fpath, [Ffile 'context.tsv']);
				case '_labeling.jpg'
					pathDest{iCount+1} = fullfile(Fpath, [Ffile 'labeling.jpg']);
				otherwise % Copy side-car name unchanged
					pathDest{iCount+1} = fullfile(Fpath, [Ffile bidsPar.sidecarName{iCar}]);
			end
            iCount = iCount+1;
        end
    end

end



