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
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

    SlicesName  = {'TraSlices' 'CorSlices' 'SagSlices'};
    for iSl=1:length(SlicesName)
        if  isfield(x.S,SlicesName{iSl})
            x.S     = rmfield(x.S,SlicesName{iSl});
        end
    end
end
