function xASL_doc_Initialize
%xASL_doc_Initialize Script to call the separate documentation crawler
% functions.
%
% FORMAT:       xASL_docu_Initialize
% 
% INPUT:        None
%
% OUTPUT:       None
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function calls individual documentation crawlers to
%               create the function description markdown files.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      xASL_doc_Initialize
% __________________________________
% Copyright 2015-2020 ExploreASL


% Create the functions markdown file
xASL_doc_Crawler('M:\SoftwareDevelopment\MATLAB\m.stritt\ExploreASL\Functions', 'M:\SoftwareDevelopment\MATLAB\m.stritt\ExploreASL\Functions\README.md')







