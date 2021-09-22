function jsonOut = xASL_bids_BIDSifyAnatJSON(jsonIn,studyPar)
%xASL_bids_BIDSifyAnatJSON Goes through the JSON structure of an anatomical file and makes sure that all the necessary conversions and checks to BIDS format are applied
%
% FORMAT: jsonOut = xASL_bids_BIDSifyAnatJSON(jsonIn)
%
% INPUT:
%   jsonIn    - JSON with the input fields (REQUIRED)
%   studyPar  - JSON with the BIDS studyPar infos (REQUIRED)
%
% OUTPUT: 
%   jsonOut - ordered and checked JSON structure
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:
% It makes all the conversions to a proper BIDS structure, checks the existence of all BIDS fields, removes superfluous fields, checks all the conditions and orderes
% the structure on the output. It works according to the normal BIDS, or ASL-BIDS definition
%
% EXAMPLE: n/a
%
% __________________________________
% Copyright 2015-2021 ExploreASL


	% Fix jsonIn based on studyPar
	jsonIn = xASL_bids_MergeStudyPar(jsonIn,studyPar,'anat');

	% Default output struct
	jsonOut = jsonIn;

	% Check fields
	if isfield(jsonIn,'RepetitionTime') && isfield(jsonIn,'RepetitionTimePreparation')

		% Make sure the size is equal
		if ~isequal(size(jsonIn.RepetitionTime),size(jsonIn.RepetitionTimePreparation))
			warning('RepetitionTime and RepetitionTimePreparation dimensions do not match...');
		else
		
			% If RepetitionTimePreparation is equal to RepetitionTime, then remove RepetitionTimePreparation
			if isnear(jsonIn.RepetitionTime,jsonIn.RepetitionTimePreparation)
				jsonOut = rmfield(jsonOut,'RepetitionTimePreparation');
			end

		end
	end
	
end
