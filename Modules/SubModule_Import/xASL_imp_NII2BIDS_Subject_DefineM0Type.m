function [jsonLocal, bJsonLocalM0isFile] = xASL_imp_NII2BIDS_Subject_DefineM0Type(studyPar, bidsPar, jsonLocal, pathM0, linkM0prefix)
%xASL_imp_NII2BIDS_Subject_DefineM0Type Define M0 Type.
%
% FORMAT: [jsonLocal, bJsonLocalM0isFile] = xASL_imp_NII2BIDS_Subject_DefineM0Type(studyPar, bidsPar, jsonLocal, pathM0, linkM0prefix)
% 
% INPUT:
%   studyPar            - JSON file with the BIDS parameters relevant for the whole study (STRUCT, REQUIRED)
%   bidsPar             - Output of xASL_imp_Config (STRUCT, REQUIRED)
%   jsonLocal           - jsonLocal struct (STRUCT, REQUIRED)
%   pathM0              - path to the M0 file (CHAR ARRAY, PATH, REQUIRED)
%   linkM0prefix        - part of the link to M0 for the intendedFor field, containing the 'perf/subject_session' (CHAR ARRAY, REQUIRED)
%
% OUTPUT:
%   jsonLocal      - jsonLocal struct (STRUCT, REQUIRED)
%   bJsonLocalM0isFile    - JSON local M0 is file (BOOLEAN)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Define M0 Type
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     [jsonLocal, bJsonLocalM0isFile] = xASL_imp_NII2BIDS_Subject_DefineM0Type(studyPar, bidsPar, jsonLocal, pathM0, linkM0prefix);
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Type of an M0 image
    bJsonLocalM0isFile = 0;
	if isfield(studyPar,'M0Type') && ~isempty(studyPar.M0Type)
		% M0Type is defined according to BIDS
		switch(studyPar.M0Type)
			case bidsPar.strM0Included
				jsonLocal.M0 = true;
			case bidsPar.strM0Separate
				jsonLocal.M0 = linkM0prefix;
				bJsonLocalM0isFile = 1;
			case bidsPar.strM0Absent
				jsonLocal.M0 = false;
			case bidsPar.strM0Estimate
				if ~isfield(studyPar,'M0Estimate') || isempty(studyPar.M0Estimate) || ~isnumeric(studyPar.M0Estimate)
					error('For M0Type Estimate, the M0Estimate numeric value has to be provided');
				end
				jsonLocal.M0Estimate = studyPar.M0Estimate;
				jsonLocal.M0 = studyPar.M0Estimate;
			otherwise
				error('Unknown value in BIDS fields M0Type');
		end
		jsonLocal.M0Type = studyPar.M0Type; % Copy the parameter
	elseif ~isfield(studyPar,'M0') || isempty(studyPar.M0) || strcmpi(studyPar.M0,'separate_scan')
		% M0-field is not defined or says 'separate_scan'
		% Look for M0 standalone or inside the sequence
		if isfield(studyPar,'M0PositionInASL4D') && (max(studyPar.M0PositionInASL4D(:))>0)
			% M0 within main ASL sequence according to parameters
            jsonLocal.M0 = true;
            jsonLocal.M0Type = bidsPar.strM0Included;
        elseif xASL_exist(pathM0)
			% Separate M0-file
            jsonLocal.M0 = linkM0prefix;
            jsonLocal.M0Type = bidsPar.strM0Separate;
            bJsonLocalM0isFile = 1;
		elseif ~isempty(strfind(jsonLocal.ASLContext,bidsPar.strM0scan))
			% M0 within the main ASL sequence according to ASL context
			jsonLocal.M0 = true;
			jsonLocal.M0Type = bidsPar.strM0Included;
		else
			% M0 completely missing
			jsonLocal.M0 = false;
			jsonLocal.M0Type = bidsPar.strM0Absent;
		end
	elseif strcmpi(studyPar.M0,'UseControlAsM0')
		% Defined to use control as M0 -> for BIDS this means M0 is absent
		jsonLocal.M0 = bidsPar.strM0Absent;
	elseif strcmpi(studyPar.M0,'no_background_suppression')
		% Defined that no Bsup is used -> for BIDS this means M0 is absent
		jsonLocal.M0 = bidsPar.strM0Absent;
	else
		jsonLocal.M0 = studyPar.M0;
		if isnumeric(studyPar.M0)
			% M0 is numeric, then define M0-estimate
			jsonLocal.M0Type = bidsPar.strM0Estimate;
			jsonLocal.M0Estimate = studyPar.M0;
		elseif xASL_exist(pathM0)
			% M0 is a path that exist - then define a separate M0-file
			jsonLocal.M0Type = bidsPar.strM0Separate;
		elseif ~isfield(jsonLocal, 'ASLContext')
			warning('jsonLocal.ASLContext missing, this may crash');
		elseif ~isempty(strfind(jsonLocal.ASLContext,bidsPar.strM0scan))
			jsonLocal.M0Type = bidsPar.strM0Included;
		else
			jsonLocal.M0Type = bidsPar.strM0Absent;
		end
	end

    
end

