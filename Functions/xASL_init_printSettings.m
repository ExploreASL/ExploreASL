function xASL_init_printSettings(x)
%xASL_init_printSettings Print chosen settings
%
% FORMAT: xASL_init_printSettings(x)
%
% INPUT:
%   x       - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:
%   n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Print chosen settings.
%
% EXAMPLE:     This is part of the initialization workflow. Check out the usage there.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:  n/a
%
% Copyright 2015-2021 ExploreASL


    %% Fallbacks
    dcm2nii = '';
    nii2bids = '';
    anonymize = '';
    bids2legacy = '';
    Structural = '';
    ASL = '';
    Population = '';
    
    %% Texts
    if x.opts.ImportModules(1)==1,   dcm2nii = 'DCM2NII ';            end
    if x.opts.ImportModules(2)==1,   nii2bids = 'NII2BIDS ';          end
    if x.opts.ImportModules(3)==1,   anonymize = 'ANONYMIZE ';        end
    if x.opts.ImportModules(4)==1,   bids2legacy = 'BIDS2LEGACY ';    end
    if x.opts.ProcessModules(1)==1,  Structural = 'Structural ';      end
    if x.opts.ProcessModules(2)==1,  ASL = 'ASL ';                    end
    if x.opts.ProcessModules(3)==1,  Population = 'Population ';      end
    
    %% Printing
    fprintf('==================================== ExploreASL Settings =====================================\n');
    if length(x.opts.DatasetRoot)>70
        fprintf('Dataset Root        ...%s\n', x.opts.DatasetRoot(end-70:end));
    else
        fprintf('Dataset Root        %s\n', x.opts.DatasetRoot); 
    end
    fprintf('Import Modules      %s%s%s%s\n', dcm2nii, nii2bids, anonymize, bids2legacy);
    fprintf('Process Modules     %s%s%s\n', Structural, ASL, Population);
    if x.opts.bPause==1
        fprintf('bPause              %s\n', 'True');
    else
        fprintf('bPause              %s\n', 'False');
    end
    fprintf('iWorker             %d\n', x.opts.iWorker);
    fprintf('nWorkers            %d\n', x.opts.nWorkers);
    fprintf('==============================================================================================\n');

end


