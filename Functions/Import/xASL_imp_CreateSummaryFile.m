function xASL_imp_CreateSummaryFile(thisSubject, PrintDICOMFields, x, bReportTotals)
%xASL_imp_CreateSummaryFile Create summary file.
%
% FORMAT: xASL_imp_CreateSummaryFile(thisSubject, PrintDICOMFields, x)
% 
% INPUT:
%   thisSubject            - Current subject struct (REQUIRED, STRUCT)
%   PrintDICOMFields       - Print DICOM fields (REQUIRED, CELL ARRAY)
%   x                      - ExploreASL x structure (REQUIRED, STRUCT)
%   x.modules.import.imPar - JSON file with structure with import parameters (REQUIRED, STRUCT)
%   bReportTotals          - true for reporting totals (OPTIONAL, DEFAULT =
%                            false)
%
% OUTPUT:
%   n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Create summary file.
%
% 1. Create summary file
% 2. Report totals (if requested)
% 3. Close the file (for writing by Matlab)
%
% For the detailed description of the overview sub-structure (thisSubject & thisVisit)
% please check out the description within xASL_imp_DetermineSubjectStructure.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     xASL_imp_CreateSummaryFile(thisSubject, PrintDICOMFields, x);
% __________________________________
% Copyright 2015-2022 ExploreASL

    
    %% Admin
    if nargin<4 || isempty(bReportTotals) || ~islogical(bReportTotals)
        bReportTotals = false;
    end

    %% 1. Create summary file
	summary_filepath = fullfile(x.modules.import.imPar.DerivativesRoot, 'ExploreASL', 'log', ['importSummary_' x.modules.import.dateTime '.tsv']);        
    summaryFileExisted = exist(summary_filepath, 'file');
    fid_summary = xASL_fOpenClose(summary_filepath, 1);

	
    %% Print headers for parameters obtained from NIfTI file
    if ~summaryFileExisted
	    fprintf(fid_summary,'subject\tvisit\tsession\tscanID\tlegacyNIfTI_filename\tdx\tdy\tdz\tdt\tnx\tny\tnz\tnt'); % \t = tab-separator

        % Print headers for parameters obtained from DICOM file
	    if exist('PrintDICOMFields','var')
		    for iField=1:length(PrintDICOMFields)
			    fprintf(fid_summary,['\t' PrintDICOMFields{iField}]);
		    end
	    end
	    fprintf(fid_summary,'\n');
    end
    

    %% Print subject
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
                        fprintf(fid_summary,'"%s"\t"%s"\t"%s"\t"%s"%s\t\n', ... % \t = tab-separator
                            x.modules.import.listsIDs.subjectIDs{iSubject}, ...
                            thisSubject.visitIDs{iVisit}, ...
                            x.modules.import.imPar.sessionNames{iSession}, ....
                            thisVisit.scanNames{iScan}, ...
                            thisSubject.summary_lines{iSubject, iVisit, iSession, iScan});
                    end
                end
            end
        end
    end
% 	fprintf(fid_summary,'\n');
	
	%% 2. Report totals
if bReportTotals
	nMissing = sum(thisSubject.globalCounts.missing_scans(:));
	nSkipped = sum(thisSubject.globalCounts.skipped_scans(:));

	% header first
	fprintf(fid_summary,'\n');
	fprintf(fid_summary,'\nSubjects\tnConverted');
	fprintf(fid_summary,['\tMissingScans (n=' num2str(nMissing) ')']);
	fprintf(fid_summary,['\tSkippedScans (n=' num2str(nSkipped) ')\n']);
    
	% then subjects row-by-row
	for iSubject=1:x.modules.import.nSubjects
		for iVisit=1:thisSubject.nVisits
            % Get fieldname
            visitFieldName = ['visit_' num2str(iVisit,'%03.f')];
            % Get visit
            thisVisit = thisSubject.(visitFieldName);
            
			fprintf(fid_summary,'"%s"', [x.modules.import.listsIDs.subjectIDs{iSubject} thisSubject.visitIDs{iVisit}]);
			fprintf(fid_summary,'\t%d',sum(thisSubject.globalCounts.converted_scans(iSubject,:,:,:)));
			
			for iSession=1:thisVisit.nSessions
				fprintf(fid_summary,'\t"');
				fprintf(fid_summary,'%s ',thisVisit.scanNames{logical(thisSubject.globalCounts.missing_scans(iSubject, iVisit, iSession,:))});
				fprintf(fid_summary,'"');
			end
			
			for iSession=1:thisVisit.nSessions
				fprintf(fid_summary,'\t"');
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
end

	%% 3. Close the file
	fclose(fid_summary);
    
    
end