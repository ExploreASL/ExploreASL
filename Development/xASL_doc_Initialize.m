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



    %% Workflow

    % Initialize ExploreASL
    x = ExploreASL_Initialize([],0);

    % Copy and modify the index README
    copyfile(fullfile(x.MyPath,'README.md'),fullfile(x.MyPath,'Development','Documentation_GitHub','index.md'));

    % Create the functions markdown file
    xASL_doc_Crawler(fullfile(x.MyPath,'Functions'), fullfile(x.MyPath,'Development','Documentation_GitHub','Functions.md'));






end
