% Compile all for Octave
% make sure we are in the src/mex subfolder of ExploreASL

MexDir    = cd;
cd(MexDir);
MexList   = xASL_adm_GetFileList(MexDir,'.*\.(c|cpp)$');
for iL=1:length(MexList)
    mex(MexList{iL});
end
