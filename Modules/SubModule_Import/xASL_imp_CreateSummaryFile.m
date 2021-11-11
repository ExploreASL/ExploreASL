function xASL_imp_CreateSummaryFile(thisSubject, imPar, PrintDICOMFields, x)
%xASL_imp_CreateSummaryFile Create summary file.
%
% FORMAT: xASL_imp_CreateSummaryFile(thisSubject, imPar, PrintDICOMFields, x)
% 
% INPUT:
%   thisSubject       - Current subject struct (REQUIRED, STRUCT)
%   imPar             - JSON file with structure with import parameters (REQUIRED, STRUCT)
%   PrintDICOMFields  - Print DICOM fields (REQUIRED, CELL ARRAY)
%   x                 - ExploreASL x structure (REQUIRED, STRUCT)
%
% OUTPUT:
%   n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Create summary file.
%
% 1. Create summary file
% 2. Report totals
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     xASL_imp_CreateSummaryFile(thisSubject, imPar, PrintDICOMFields, x);
% __________________________________
% Copyright 2015-2021 ExploreASL

    
    %% 1. Create summary file
	summary_filepath = fullfile(imPar.TempRoot, ['import_summary_' thisSubject.name '.csv']);
	fid_summary = fopen(summary_filepath,'wt');
	
    % Print headers for parameters obtained from NIfTI file
	fprintf(fid_summary,'subject,visit,session,scan,filename,dx,dy,dz,dt,nx,ny,nz,nt');
	
    % Print headers for parameters obtained from DICOM file
	if exist('PrintDICOMFields','var')
		for iField=1:length(PrintDICOMFields)
			fprintf(fid_summary,[',' PrintDICOMFields{iField}]);
		end
	end
	fprintf(fid_summary,'\n');
    
    for iSubject=1:x.modules.import.nSubjects
        for iVisit=1:thisSubject.nVisits
            % Get fieldname
            visitFieldName = ['visit_' num2str(iVisit,'%03.f')];
            % Get visit
            thisVisit = thisSubject.(visitFieldName);
            for iScan=1:thisVisit.nScans
                for iSession=1:thisVisit.nSessions
                    if thisSubject.globalCounts.converted_scans(iSubject, iVisit, iSession, iScan) || ...
                            thisSubject.globalCounts.skipped_scans(iSubject, iVisit, iSession, iScan) || ...
                            thisSubject.globalCounts.missing_scans(iSubject, iVisit, iSession, iScan)
                        fprintf(fid_summary,'"%s","%s","%s","%s"%s,\n', ...
                            x.modules.import.listsIDs.subjectIDs{iSubject}, ...
                            thisSubject.visitIDs{iVisit}, ...
                            imPar.sessionNames{iSession}, ....
                            thisVisit.scanNames{iScan}, ...
                            thisSubject.summary_lines{iSubject, iVisit, iSession, iScan});
                    end
                end
            end
        end
    end
	fprintf(fid_summary,'\n');
	
	nMissing = sum(thisSubject.globalCounts.missing_scans(:));
	nSkipped = sum(thisSubject.globalCounts.skipped_scans(:));
	
	%% 2. Report totals
    
	% header first
	fprintf(fid_summary,'\n');
	fprintf(fid_summary,'\nSubject,nConverted');
	fprintf(fid_summary,[',MissingScans (n=' num2str(nMissing) ')']);
	fprintf(fid_summary,[',SkippedScans (n=' num2str(nSkipped) ')\n']);
    
	% then subjects row-by-row
	for iSubject=1:x.modules.import.nSubjects
		for iVisit=1:thisSubject.nVisits
            % Get fieldname
            visitFieldName = ['visit_' num2str(iVisit,'%03.f')];
            % Get visit
            thisVisit = thisSubject.(visitFieldName);
            
			fprintf(fid_summary,'"%s"', [x.modules.import.listsIDs.subjectIDs{iSubject} thisSubject.visitIDs{iVisit}]);
			fprintf(fid_summary,',%d',sum(thisSubject.globalCounts.converted_scans(iSubject,:,:,:)));
			
			for iSession=1:thisVisit.nSessions
				fprintf(fid_summary,',"');
				fprintf(fid_summary,'%s ',thisVisit.scanNames{logical(thisSubject.globalCounts.missing_scans(iSubject, iVisit, iSession,:))});
				fprintf(fid_summary,'"');
			end
			
			for iSession=1:thisVisit.nSessions
				fprintf(fid_summary,',"');
				fprintf(fid_summary,'%s ',thisVisit.scanNames{logical(thisSubject.globalCounts.skipped_scans(iSubject, iVisit, iSession,:))});
				fprintf(fid_summary,'"');
			end
			
			fprintf(fid_summary,'\n');
		end
	end
	
	% and a grand total of missing and skipped
	if nMissing>0
		fprintf(2,'Number of missing scans: %d\n',nMissing);
	end
	if nSkipped>0
		fprintf(2,'Number of skipped scans: %d\n',nSkipped);
	end
	fclose(fid_summary);
    
    
end

