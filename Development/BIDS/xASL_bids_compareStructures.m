function [identical,results] = xASL_bids_compareStructures(rootStructureA,rootStructureB)
%xASL_bids_compareStructures Function that compares two BIDS folders with several subfolders and studies and prints the differences.
%
% FORMAT: results = xASL_bids_compareStructures(rootStructureA,rootStructureB)
%
% INPUT:
%        rootStructureA     - path to first BIDS structure (REQUIRED)
%        rootStructureB     - path to second BIDS structure (REQUIRED)
%
% OUTPUT:
%        identical          - Returns 1 if both folder structures are identical and 0 if not
%        results            - structure containing (possible) differences of both folder structures
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      Function that compares two BIDS folders with several subfolders and studies and prints the differences.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          rootStructureA = '...\bids-examples\ds001';
%                   rootStructureB = '...\bids-examples\ds001_exact_copy'
%                   [identical,results] = xASL_bids_compareStructures(rootStructureA,rootStructureB);
%
% REFERENCES:       ...
% __________________________________
% Copyright @ 2015-2020 ExploreASL




