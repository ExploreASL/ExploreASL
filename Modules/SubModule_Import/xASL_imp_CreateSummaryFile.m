function xASL_imp_CreateSummaryFile(imPar, numOf, listsIDs, PrintDICOMFields, converted_scans, skipped_scans, missing_scans, scanNames, summary_lines, fid_summary)
%xASL_imp_CreateSummaryFile Create summary file.
%
% FORMAT: xASL_imp_CreateSummaryFile(imPar, PrintDICOMFields, converted_scans, skipped_scans, missing_scans, subjectIDs, visitIDs, scanNames, summary_lines, nSubjects, nVisits, nSessions, fid_summary)
% 
% INPUT:
%   imPar             - JSON file with structure with import parameters (REQUIRED, STRUCT)
%   numOf             - Struct defining the number of subjects, visists, sessions, scans etc. (REQUIRED, STRUCT)
%   listsIDs          - Struct defining different ID lists (REQUIRED, STRUCT)
%   PrintDICOMFields  - Print DICOM fields (REQUIRED, CELL ARRAY)
%   converted_scans   - Converted scans (REQUIRED)
%   skipped_scans     - Skipped scans (REQUIRED)
%   missing_scans     - Missing scans (REQUIRED)
%   scanNames         - Scan names (REQUIRED, CELL ARRAY)
%   summary_lines     - Summary lines (REQUIRED, CELL ARRAY)
%   fid_summary       - File ID summary (REQUIRED, INTEGER)
%
% OUTPUT:
%   n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Create summary file.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     xASL_imp_CreateSummaryFile(imPar, PrintDICOMFields, converted_scans, skipped_scans, missing_scans, subjectIDs, visitIDs, scanNames, summary_lines, nSubjects, nVisits, nSessions, fid_summary);
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Extract structs
    
    % Extract ID struct
    subjectIDs = listsIDs.subjectIDs;
    visitIDs = listsIDs.visitIDs;
    
    % Extract number of struct
    nSubjects = numOf.nSubjects;
    nVisits = numOf.nVisits;
    nSessions = numOf.nSessions;
    nScans = numOf.nScans;

    
    %% Create summary file
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
	
	for iScan=1:nScans
		for iSubject=1:nSubjects
			for iVisit=1:nVisits
				for iSession=1:nSessions
					if converted_scans(iSubject, iVisit, iSession, iScan) || skipped_scans(iSubject, iVisit, iSession, iScan) || missing_scans(iSubject, iVisit, iSession, iScan)
						fprintf(fid_summary,'"%s","%s","%s","%s"%s,\n', subjectIDs{iSubject}, visitIDs{iVisit}, imPar.sessionNames{iSession}, scanNames{iScan}, summary_lines{iSubject, iVisit, iSession, iScan});
					end
				end
			end
		end
	end
	fprintf(fid_summary,'\n');
	
	nMissing = sum(missing_scans(:));
	nSkipped = sum(skipped_scans(:));
	
	%% Report totals
    
	% header first
	fprintf(fid_summary,'\n');
	fprintf(fid_summary,'\nSubject,nConverted');
	fprintf(fid_summary,[',MissingScans (n=' num2str(nMissing) ')']);
	fprintf(fid_summary,[',SkippedScans (n=' num2str(nSkipped) ')\n']);
    
	% then subjects row-by-row
	for iSubject=1:nSubjects
		for iVisit=1:nVisits
			fprintf(fid_summary,'"%s"', [subjectIDs{iSubject} visitIDs{iVisit}]);
			fprintf(fid_summary,',%d',sum(converted_scans(iSubject,:,:,:)));
			
			for iSession=1:nSessions
				fprintf(fid_summary,',"');
				fprintf(fid_summary,'%s ',scanNames{logical(missing_scans(iSubject, iVisit, iSession,:))});
				fprintf(fid_summary,'"');
			end
			
			for iSession=1:nSessions
				fprintf(fid_summary,',"');
				fprintf(fid_summary,'%s ',scanNames{logical(skipped_scans(iSubject, iVisit, iSession,:))});
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



    
    