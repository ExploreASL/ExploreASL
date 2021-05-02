function xASL_imp_CreateSummaryFile(imPar, PrintDICOMFields, x, scanNames, summary_lines, fid_summary)
%xASL_imp_CreateSummaryFile Create summary file.
%
% FORMAT: xASL_imp_CreateSummaryFile(imPar, numOf, listsIDs, PrintDICOMFields, globalCounts, scanNames, summary_lines, fid_summary)
% 
% INPUT:
%   imPar             - JSON file with structure with import parameters (REQUIRED, STRUCT)
%   PrintDICOMFields  - Print DICOM fields (REQUIRED, CELL ARRAY)
%   globalCounts      - Converted, skipped & missing scans (REQUIRED, STRUCT)
%   scanNames         - Scan names (REQUIRED, CELL ARRAY)
%   summary_lines     - Summary lines (REQUIRED, CELL ARRAY)
%   fid_summary       - File ID summary (REQUIRED, INTEGER)
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
% EXAMPLE:     xASL_imp_CreateSummaryFile(imPar, numOf, listsIDs, PrintDICOMFields, globalCounts, scanNames, summary_lines, fid_summary);
% __________________________________
% Copyright 2015-2021 ExploreASL

    
    %% 1. Create summary file
    
    numOf = x.modules.import.numOf;
    listsIDs = x.modules.import.listsIDs;
    globalCounts = x.modules.import.globalCounts;
    
	summary_filepath = fullfile(imPar.AnalysisRoot, 'import_summary.csv');
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
	
	for iScan=1:numOf.nScans
		for iSubject=1:numOf.nSubjects
			for iVisit=1:numOf.nVisits
				for iSession=1:numOf.nSessions
					if globalCounts.converted_scans(iSubject, iVisit, iSession, iScan) || globalCounts.skipped_scans(iSubject, iVisit, iSession, iScan) || globalCounts.missing_scans(iSubject, iVisit, iSession, iScan)
						fprintf(fid_summary,'"%s","%s","%s","%s"%s,\n', listsIDs.subjectIDs{iSubject}, listsIDs.visitIDs{iVisit}, imPar.sessionNames{iSession}, scanNames{iScan}, summary_lines{iSubject, iVisit, iSession, iScan});
					end
				end
			end
		end
	end
	fprintf(fid_summary,'\n');
	
	nMissing = sum(globalCounts.missing_scans(:));
	nSkipped = sum(globalCounts.skipped_scans(:));
	
	%% 2. Report totals
    
	% header first
	fprintf(fid_summary,'\n');
	fprintf(fid_summary,'\nSubject,nConverted');
	fprintf(fid_summary,[',MissingScans (n=' num2str(nMissing) ')']);
	fprintf(fid_summary,[',SkippedScans (n=' num2str(nSkipped) ')\n']);
    
	% then subjects row-by-row
	for iSubject=1:numOf.nSubjects
		for iVisit=1:numOf.nVisits
			fprintf(fid_summary,'"%s"', [listsIDs.subjectIDs{iSubject} listsIDs.visitIDs{iVisit}]);
			fprintf(fid_summary,',%d',sum(globalCounts.converted_scans(iSubject,:,:,:)));
			
			for iSession=1:numOf.nSessions
				fprintf(fid_summary,',"');
				fprintf(fid_summary,'%s ',scanNames{logical(globalCounts.missing_scans(iSubject, iVisit, iSession,:))});
				fprintf(fid_summary,'"');
			end
			
			for iSession=1:numOf.nSessions
				fprintf(fid_summary,',"');
				fprintf(fid_summary,'%s ',scanNames{logical(globalCounts.skipped_scans(iSubject, iVisit, iSession,:))});
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



    
    