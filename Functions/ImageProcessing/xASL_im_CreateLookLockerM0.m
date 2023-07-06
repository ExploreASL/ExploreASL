function x = xASL_im_CreateLookLockerM0(x)
%xASL_im_CreateLookLockerM0 Create M0 out of Look-Locker control images using FSL asl_calib
%
% FORMAT: x = xASL_im_CreateLookLockerM0(x)
%
% INPUT:
%   x           - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
% OUTPUT:
%   x           - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%                 assigns the following values:
%                 x.Path_M0_json - Path to newly create M0 json
%-----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function creates an M0 out of Look-Locker control images using FSL asl_calib

PathControl4D = fullfile(x.dir.SESSIONDIR ,'Control4D_FSLInput.nii');
asl_calib_directory = fullfile(x.dir.SESSIONDIR ,'asl_calib_output/');
asl_calib_M0 = fullfile(asl_calib_directory ,'M0t.nii');
if ~exist(asl_calib_directory,'dir')
    mkdir(asl_calib_directory)
end
% add LL recovery correction !!!!!!!!!!!!!!
TIstring = regexprep(num2str((x.Q.Initial_PLD(1:end/2)/1000)'),' +', ',');
ASL_Calib_Options = ['-c ' PathControl4D ' --mode satrecov --tissref wm --te  13.894  -m ' x.P.Path_PVwm ' --tis ' TIstring ' --fa ' num2str(x.Q.FlipAngle) ' -o ' asl_calib_directory];
[~, resultFSL] = xASL_fsl_RunFSL(['asl_calib ' ASL_Calib_Options], x);
xASL_Copy(asl_calib_M0,x.P.Path_M0);

% create M0.json
xASL_Copy(x.P.Path_ASL4D_json,x.P.Path_M0_json);
NewM0Json = xASL_io_ReadJson(x.P.Path_M0_json);
NewM0Json.FlipAngle = 90; % change for newly created M0
NewM0Json.PostLabelingDelay = []; % remove
NewM0Json.RepetitionTimePreparation = 3; % !!! CHANGE !!!! This is due to an error in import %
xASL_io_WriteJson(x.P.Path_M0_json,NewM0Json,1);
end