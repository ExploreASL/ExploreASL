function xASL_imp_CreateSummaryFile(thisSubject, PrintDICOMFields, x)
%xASL_imp_CreateSummaryFile Create summary file.
%
% FORMAT: xASL_imp_CreateSummaryFile(thisSubject, PrintDICOMFields, x)
% 
% INPUT:
%   thisSubject            - Current subject struct (REQUIRED, STRUCT)
%   PrintDICOMFields       - Print DICOM fields (REQUIRED, CELL ARRAY)
%   x                      - ExploreASL x structure (REQUIRED, STRUCT)
%   x.modules.import.imPar - JSON file with structure with import parameters (REQUIRED, STRUCT)
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
% For the detailed description of the overview sub-structure (thisSubject & thisVisit)
% please check out the description within xASL_imp_DetermineSubjectStructure.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     xASL_imp_CreateSummaryFile(thisSubject, PrintDICOMFields, x);
% __________________________________
% Copyright 2015-2022 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________


    
    %% 1. Create summary file
	summary_filepath = fullfile(x.modules.import.imPar.DerivativesRoot, 'ExploreASL', 'log', ['import_summary_sub-' thisSubject.name '.csv']);
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
    
    % Check if the converted_scans, skipped_scans, and missing_scans tables
    % have equal sizes, otherwise try to fix it

	% We are expecting a 4D matrix. So we have to ask for 4 dimensions as this is required on the lines below
	% For example size(ones(5,1,1,1)      = [5 1]
	% but         size(ones(5,1,1,1),1:4) = [5 1 1 1]
    sizeConvertedScans = size(thisSubject.globalCounts.converted_scans, 1:4);
    sizeSkippedScans = size(thisSubject.globalCounts.skipped_scans, 1:4);
    sizeMissingScans = size(thisSubject.globalCounts.missing_scans, 1:4);

    if ~isequal(sizeConvertedScans, sizeSkippedScans)
        warning('Skipped scans had a different size than converted scans, fixing this but summary file could be incorrect');
        skippedScans = uint8(zeros(size(thisSubject.globalCounts.converted_scans))); % by default assume that scans are not skipped (zeros), unless this was set
        skippedScans(logical(thisSubject.globalCounts.skipped_scans)) = 1;
        thisSubject.globalCounts.skipped_scans = skippedScans;
    end
    if ~isequal(sizeConvertedScans, sizeMissingScans)
        warning('Missing scans had a different size than converted scans, fixing this but summary file could be incorrect');
        missingScans = uint8(ones(size(thisSubject.globalCounts.converted_scans))); % by default we assume that scans are missing (ones), unless we confirmed their presence
        missingScans(~logical(thisSubject.globalCounts.missing_scans)) = 0;
        thisSubject.globalCounts.missing_scans = missingScans;
    end

    for iSubject=1:x.modules.import.nSubjects
        for iVisit=1:thisSubject.nVisits
            % Get fieldname
            visitFieldName = ['visit_' num2str(iVisit,'%03.f')];
            % Get visit
            thisVisit = thisSubject.(visitFieldName);
            for iScan=1:thisVisit.nScans
                for iSession=1:thisVisit.nSessions

                    % Here we try to skip this table if something went wrong in dicom2nii
                    bSkipIt = false;
                    if x.modules.import.nSubjects>sizeConvertedScans(1)
                        warning('Something went wrong with number of subjects');
                        bSkipIt = true;
                    elseif thisSubject.nVisits>sizeConvertedScans(2)
                        warning('Something went wrong with number of visits');
                        bSkipIt = true;
                    elseif thisVisit.nSessions>sizeConvertedScans(3)
                        warning('Something went wrong with number of sessions');
                        bSkipIt = true;
                    elseif thisVisit.nScans>sizeConvertedScans(4)
                        warning('Something went wrong with number of scans');
                        bSkipIt = true;
                    end


                    if ~bSkipIt && ( thisSubject.globalCounts.converted_scans(iSubject, iVisit, iSession, iScan) || ...
                            thisSubject.globalCounts.skipped_scans(iSubject, iVisit, iSession, iScan) || ...
                            thisSubject.globalCounts.missing_scans(iSubject, iVisit, iSession, iScan) )
                        fprintf(fid_summary,'"%s","%s","%s","%s"%s,\n', ...
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

