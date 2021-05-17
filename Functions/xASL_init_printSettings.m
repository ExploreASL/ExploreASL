%% -----------------------------------------------------------------------
%% Print chosen settings
function xASL_init_printSettings(x)

    % Fallbacks
    dcm2nii = '';
    nii2bids = '';
    anonymize = '';
    bids2legacy = '';
    Structural = '';
    ASL = '';
    Population = '';
    
    % Texts
    if x.ImportModules(1)==1,   dcm2nii = 'DCM2NII';            end
    if x.ImportModules(2)==1,   nii2bids = 'NII2BIDS';          end
    if x.ImportModules(3)==1,   anonymize = 'ANONYMIZE';        end
    if x.ImportModules(4)==1,   bids2legacy = 'BIDS2LEGACY';    end
    if x.ProcessModules(1)==1,  Structural = 'Structural';      end
    if x.ProcessModules(2)==1,  ASL = 'ASL';                    end
    if x.ProcessModules(3)==1,  Population = 'Population';      end
    
    % Printing
    fprintf('==================================== ExploreASL Settings =====================================\n');
    if length(x.DataParPath)>70,    fprintf('DataParPath         ...%s\n', x.DataParPath(end-70:end));
    else,                           fprintf('DataParPath         %s\n', x.DataParPath); end
    fprintf('Import Modules      %s %s %s %s\n', dcm2nii, nii2bids, anonymize, bids2legacy);
    fprintf('Process Modules     %s %s %s\n', Structural, ASL, Population);
    if x.bPause==1
        fprintf('bPause              %s\n', 'True');
    else
        fprintf('bPause              %s\n', 'False');
    end
    fprintf('iWorker             %d\n', x.iWorker);
    fprintf('nWorkers            %d\n', x.nWorkers);
    fprintf('==============================================================================================\n');

end


