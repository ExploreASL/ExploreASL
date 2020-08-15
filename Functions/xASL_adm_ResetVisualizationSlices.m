function [x] = xASL_adm_ResetVisualizationSlices(x)
%xASL_adm_ResetVisualizationSlices Removes any predefined slices that should be visualized, 
% allowing to show the default slices. Comes in handy when different
% pipeline visualization parts are repeated.
%
% FORMAT:       [x] = xASL_adm_ResetVisualizationSlices(x)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Removes any predefined slices that should be visualized, 
%               allowing to show the default slices. Comes in handy when different
%               pipeline visualization parts are repeated.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL

    SlicesName  = {'TraSlices' 'CorSlices' 'SagSlices'};
    for iSl=1:length(SlicesName)
        if  isfield(x.S,SlicesName{iSl})
            x.S     = rmfield(x.S,SlicesName{iSl});
        end
    end
end

