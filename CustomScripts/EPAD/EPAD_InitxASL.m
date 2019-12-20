function [x] = EPAD_InitxASL
%EPAD_INITXASL This function calls ExploreASL only for initialization of
%paths etc

%% START THIS FUNCTION IN THE CustomScripts/EPAD folder

CurrentDir = pwd;

if ~isempty(which(mfilename)) % current function
    addpath(fileparts(which(mfilename))); % add to path
end

% Initialize ExploreASL
if ~isempty(which('ExploreASL_Master'))
    DirExploreASL = fileparts(which('ExploreASL_Master'));
    cd(DirExploreASL);
else
    % try going two folders lower
    Fpath = fileparts(fileparts(fileparts(which(mfilename))));
    cd(Fpath);
end
    

close all;
fclose all;
clc;
x = ExploreASL_Master('',0);

cd(CurrentDir);

end

