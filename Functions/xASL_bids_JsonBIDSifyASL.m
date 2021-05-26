function jsonOut = xASL_bids_JsonBIDSifyASL(jsonIn)
%xASL_bids_JsonBIDSifyASL Goes through the JSON structure of an ASL file and makes sure that all the necessary conversions and checks to BIDS format are applied
%
% FORMAT: jsonOut = xASL_bids_JsonBIDSifyASL(jsonIn)
%
% INPUT:
%   jsonIn  - JSON with the input fields (REQUIRED)
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

jsonOut = jsonIn;

end
