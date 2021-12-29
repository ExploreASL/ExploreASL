function xASL_quant_FEAST(x)
%xASL_quant_FEAST Computation FEAST-based transit times
%
% FORMAT: xASL_quant_FEAST(x)
% 
% INPUT:
%   x           - struct containing pipeline environment parameters
% OUTPUT:
%   x           - struct containing pipeline environment parameters
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function quantifies ATT using the FEAST equations,
% using crushed and non-crushed sessions, of which the ratio is
% proportional to ATT.
% Note that the order of sessions should be 1) crushed 2) non-crushed
%
% This function runs the following steps:
%
% 1. Skip this function if no FEAST data available
% 2. Admin
% 3. Load data & correct for timing differences (PLD etc)
% 4. Smooth and clip CBF maps & FEAST ratio
% 5. Compute TT maps
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_quant_FEAST(x);
% REFERENCE: JJ Wang, 2003 MRM; Y Chen, 2012 MAGMA
% __________________________________
% Copyright 2015-2021 ExploreASL


%% -------------------------------------------------------------
%% 1. Skip this function if no FEAST data available
if ~(x.dataset.nSessions>1 && isfield(x,'session') && isfield(x.session,'options') && strcmp(x.session.options{1},'crushed') && strcmp(x.session.options{2},'non-crushed'))
    return;
elseif ~(strcmp(x.dir.SESSIONDIR(length(x.dir.SUBJECTDIR)+2:end),'ASL_2')) % Computation is performed if CurrentSession=session 2
    return;
end

%% -------------------------------------------------------------
%% 2. Admin
for iSession=1:2
    PathCBF{iSession} = fullfile(x.D.PopDir, ['qCBF_' x.P.SubjectID '_' x.SESSIONS{iSession} '.nii']);
end
if ~xASL_exist(PathCBF{1},'file') || ~xASL_exist(PathCBF{2},'file')
    return; % skip if files dont exist (yet)
end

fprintf('%s\n','Saving TT nifti');
xASL_adm_CreateDir(x.D.TTCheckDir);

%% -------------------------------------------------------------
%% 3. Load data & correct for timing differences (PLD etc)
for iSession=1:2
    % Load data
    CBF{iSession} = xASL_io_Nifti2Im(PathCBF{iSession});
    SliceNumber = xASL_io_Nifti2Im(fullfile(x.D.PopDir, ['SliceGradient_extrapolated_' x.P.SubjectID '_' x.SESSIONS{iSession} '.nii']));

	% Obtain the correct SliceReadoutTime
	SliceReadoutTime = xASL_quant_SliceTiming(x,x.P.Path_ASL4D);
	
    % Correct different PLD scales
	SliceNumber = round(SliceNumber);
	SliceNumber(SliceNumber<1) = 1;
	SliceNumber(SliceNumber>length(SliceReadoutTime)) = length(SliceReadoutTime);
	PLD{iSession} = x.Q.Initial_PLD + SliceReadoutTime(SliceNumber);
	
    CBF{iSession} = CBF{iSession}./(exp(PLD{iSession}./x.Q.BloodT1) / (2.*x.Q.LabelingEfficiency.*x.Q.BloodT1 .* (1- exp(-x.Q.LabelingDuration./x.Q.BloodT1)) ));
end
% Average different PLD scales
PLD_combined= (PLD{1}+PLD{2})./2;

%% -------------------------------------------------------------
%% 4. Smooth and clip CBF maps & FEAST ratio
for iSession=1:2
    % Masking
    CBF{iSession}(~x.S.masks.skull) = NaN; % this is not masked for WM or GM, we could do that
    % smooth maps, ignoring NaNs
    % CAVE: NaNs are interpolated with data, hence the maps should be masked later!

    CBF{iSession} = xASL_im_ndnanfilter(CBF{iSession},'gauss',[8 8 8],0);
end

FEAST_ratio = CBF{1}./CBF{2}; % crushed/non-crushed
FEAST_ratio(FEAST_ratio>1) = 1; % clip @ 1 (== ATT = PLD)
FEAST_ratio(FEAST_ratio<0) = 0; % clip @ 0 (== infinite ATT)

%% -------------------------------------------------------------
%% 5. Compute TT maps
qnt_PLDdecay = exp(-PLD_combined/x.Q.BloodT1);
qnt_combidecay = exp( (-x.Q.LabelingDuration - PLD_combined) / x.Q.BloodT1);
TT = -x.Q.BloodT1 .* reallog( FEAST_ratio .* (qnt_PLDdecay  - qnt_combidecay ) + qnt_combidecay );

xASL_io_SaveNifti(PathCBF{1}, fullfile(x.D.PopDir, ['TT_'  x.P.SubjectID '.nii']), TT, 32);


end

