function [Parms] = xASL_adm_uiGetInput(Parms)
%xASL_adm_uiGetInput Checks whether input fields are present,
% or requests them.
%
% FORMAT:       [Parms] = xASL_adm_uiGetInput(Parms)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Checks whether input fields are present, or requests them.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL

% Check if we are in display or CLI mode
 if ~usejava('desktop') || ~usejava('jvm') || ~feature('ShowFigureWindows')
     UseGUI = false;
 else
     UseGUI = true;
 end


% First check for valid input
% First check for valid input
if  isfield(Parms,'Input2Check')
    for iI=1:length(Parms.Input2Check)
        if ~isfield(Parms,Parms.Input2Check{iI})
            ExistField(iI)  = 0;
        elseif isempty(Parms.(Parms.Input2Check{iI}))
            ExistField(iI)  = 0;
        else
            ExistField(iI)  = 1;
        end
    end

    for iI=1:length(Parms.Input2Check)
        if ~ExistField(iI)
            if strcmp(Parms.InputFormat{iI},'string')
                    if UseGUI
                        TempStr = inputdlg({['Enter ' Parms.InputNaming{iI}]},Parms.Input2Check{iI});
                    else
                        TempStr = input(['Please enter ' Parms.InputNaming{iI}]);
                    end
                    if  isempty(TempStr); return; end
                    Parms.(Parms.Input2Check{iI}) = TempStr;
            elseif  strcmp(Parms.InputFormat{iI},'dir')
                    if UseGUI
                        TempDir = uigetdir([],['Select folder containing ' Parms.InputNaming{iI}]);
                    else
                        TempDir = input(['Please enter folder containing ' Parms.InputNaming{iI}]);
                    end
                    if sum(TempDir==0); return; end
                    Parms.(Parms.Input2Check{iI}) = TempDir;
            elseif  strcmp(Parms.InputFormat{iI},'file')
                    if UseGUI
                        [name, pathstr] = uigetfile('*.*',['Select ' Parms.InputNaming{iI}]);
                    else
                        pathstr = input(['Please enter ' Parms.InputNaming{iI}]);
                    end
                    if sum(pathstr==0); return; end
                    Parms.(Parms.Input2Check{iI}) = fullfile(pathstr,name);
            end
        end
    end
end






end