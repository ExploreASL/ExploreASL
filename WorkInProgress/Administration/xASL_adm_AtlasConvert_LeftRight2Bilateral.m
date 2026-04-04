function xASL_adm_AtlasConvert_LeftRight2Bilateral(pathNiftii, pathTSVin, pathTSVout, atlasType)

% Load atlas saved in pathNiftii with labels in pathTSV and convert to bilateral using the scheme specified in atlasType

% Load image and label information
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________



IM = xASL_io_Nifti2Im(pathNiftii);
TSV = xASL_tsvRead(pathTSVin);

checkNR = unique(IM(:))

switch(atlasType)

	case 'AAL3v1'
		% The left-right ROIs are averaged as ExploreASL divides ROIs in left-right-bilateral automatically,
		% except for ROIs 113-120 (vermis) and 169-170 (raphe medial and dorsal) which are bilateral already).
		% 1:112 -> 1:56
		% 113:120 -> 57:64
		% 121:168 -> 65:88
		% 169:170 -> 89:90
		TSV{1} = '1 Precentral_L 1';% Fix the first entry
		for i = 1:length(TSV)
			token = regexp(TSV{i}, '^.\d* (.*) \d*$', 'tokens');
			TSV{i} = token{1}{1};
		end
		

		for iROI = 1:56
			IM(IM==(iROI*2)-1) = iROI;
			IM(IM==iROI*2) = iROI;
			TSV{iROI} = TSV{iROI*2}(1:end-2);
		end
		for iROI = 113:120
			IM(IM==iROI) = iROI-56;
			TSV{iROI-56} = TSV{iROI};
		end
		for iROI = (121:144)-120
			IM(IM==(iROI*2)-1+120) = iROI+64;
			IM(IM==iROI*2+120) = iROI+64;
			TSV{iROI+64} = TSV{iROI*2+120}(1:end-2);
		end
		for iROI = 169:170
			IM(IM==iROI) = iROI-80;
			TSV{iROI-80} = TSV{iROI};
		end

		TSV = TSV(1:90);
	case 'WMPM_Type_III'
		%% The same for WMPM
		IMnew = zeros(size(IM));

		% Remove double _
		TSV = strrep(TSV, '__', '_');

		% Get left or right
		LeftIs = contains(TSV(:,1), 'left');
		RightIs = contains(TSV(:,1), 'right');
		TSV(LeftIs,2) = {1}; % left
		TSV(RightIs,2) = {2}; % right

		% Remove trailing left/right
		TSV(:,1) = strrep(TSV(:,1), '_left', '');
		TSV(:,1) = strrep(TSV(:,1), 'left', '');
		TSV(:,1) = strrep(TSV(:,1), '_right', '');
		TSV(:,1) = strrep(TSV(:,1), 'right', '');

		% Bugfix
		TSV{108,1} = strrep(TSV{108,1}, '_(a_part_of_MCP)', '');

		iROInew = 1;
		reportSingleROI = {''};
		for iROI=1:length(TSV)
			% is there a ROI with a similar name
			if ~isempty(TSV{iROI,1})
				indicesAre = find(strcmp(TSV(:,1), TSV{iROI,1}));

				if numel(indicesAre)==1
					reportSingleROI{end+1,1} = TSV{iROI,1};
					reportSingleROI{end,2} = iROI;
				end

				% assuming they are always either left or right
				TSVnew{iROInew} = TSV{iROI,1};
				for iIndices=1:length(indicesAre)
					IMnew(IM==indicesAre(iIndices)) = iROInew;
				end
				% Remove current ROI
				TSV(indicesAre,1) = {''};
				iROInew = iROInew+1;
			end
		end
		TSV = TSVnew';
end

% Save the image
xASL_io_SaveNifti(pathNiftii, pathNiftii, IM, []);
xASL_tsvWrite(TSV, pathTSVout, 1);
IM = uint8(IM);
save([pathNiftii '.mat'],'IM');
xASL_adm_GzipNifti(pathNiftii);

end
