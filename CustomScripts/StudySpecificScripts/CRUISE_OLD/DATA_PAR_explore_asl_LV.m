function x = DATA_PAR(x)
%DATA_PAR Part of ASL pipeline, loading basic study settings
% AMC ASL pipeline, HJMM Mutsaerts 2014

x.MYPATH             = 'L:\basic\divi\Users\lvaclavu\MATLAB\ASL_pipeline_old'; %'C:\Users\lvaclavu\MATLAB\ASL_pipeline_old';%'C:\Users\lvaclavu\MATLAB\ExploreASL'; % define the root location of the ASL toolbox here
x.DIPpath            = 'L:\basic\divi\Users\lvaclavu\MATLAB\DIPimage 2.8.1';
x.SPMpath            = 'L:\basic\divi\Users\lvaclavu\MATLAB\spm12\6225_20141001';
addpath(genpath(fullfile(x.MYPATH,'spmwrapperlib')));

% Define study
x.name               = 'CRUISE'; % name of the study
x.D.ROOT               = 'L:\basic\divi\Projects\cruise\ASL\CVR_CRUISE_T1_LE'; % root directory of study analysis data, as imported by pipeline import script
x.subject_regexp     = 'cruise102';% '^cruise0\d{2}$'%'cruise001'%'^cruise\d{3}$';
x.exclusion 		   = ''; % define exclusions, which subjects not to include in the pipeline
nSubjects                  = length(xASL_adm_GetFsList( x.D.ROOT, x.subject_regexp,1)) - length(x.exclusion);
x.SESSIONS           = {'ASL_1','ASL_2'}; %{'ASL_1'}; {'ASL_1','ASL_2'} % Define sessions, if FEAST: 1=no crushed, 2=crushed
nSessions 				   = length(x.SESSIONS);
x.session.code 	   = repmat([1:1:nSessions]',nSubjects,1);
x.session.options    = {'baseline' 'acetazolamide'}; % this is how the sessions will be called. If FEAST, this should be {'non-crushed' 'crushed'}

% current dartel folder is dual compartment, T1a 1800, LE 0.85, att 700
x.QUANT_model = 'dual_compartment' % (Wang 2002)
%x.Q.BloodT1     = 1898;						
%x.Q.LabelingEfficiency = 0.86; % remember this gets further reduced by 83% due to 2 background suppression pulses (x.Q.LabelingEfficiency*0.83)
%x.Q.ATT     = 543;%post=509; % PLUS 500 ms (Take average of Chen Wang Detre)-Average of crushed minus average of noncrushed and add the result to here
% https://link.springer.com/content/pdf/10.1007%2Fs10334-011-0276-5.pdf		
		
% Define groups here
% Here you can define variables. Name of the variable should be the same
% as the FileName
% x.group{1}.name 					  = 'ControlPatient'; % name of the group
% load( fullfile(x.D.ROOT, 'ControlPatient.mat') ); % code file of group
% x.group{1}.code 					  = ControlPatient; % load code file
% x.group{1}.options 			      = {'control' 'patient'}; % names of options (e.g. cohorts)
% x.group{1}.sets_1_2_sample  	      = 1; % for t-test, is two-sampled here, independent variables. for sessions, this is always one-sampled
% clear ControlPatient
 
% % Specify T1blood values for all subjects
    %load(fullfile(x.D.ROOT,'qnt_T1a.mat'));
    x.Q.BloodT1         = 'qnt_T1a';
% % Specify labelling efficiency values for all subjects
    %load(fullfile(x.D.ROOT,'qnt_lab_eff.mat'));
    x.Q.LabelingEfficiency     = 'qnt_lab_eff';
% % Specify ATT values for all subjects
    %load(fullfile(x.D.ROOT,'qnt_ATT.mat'));
    x.Q.ATT         = 'qnt_ATT'; 
    
% % Specify T1t values for all subjects
    %load(fullfile(x.D.ROOT,'qnt_T1t.mat'));
    %x.Q.ATT          = 'qnt_T1t';   
  
  
% list parameters here
x.M0_conventionalProcessing   = 0 % 1 uses conventional M0 processing, 0 new image processing (improved masking & smoothing). default = 0)
x.M0                          = 'separate_scan'; % no_background_suppression 
% 'separate_scan'; % no_background_suppression
% 3.7394*10^6; Intera AMC Dennis PET-ASL study 
% Ingenia AMC 
% based on AgeIV ratio: 3.7394*10^6/12.264 = 3.0491e+005
% based on Novice ASL protocol M0 measurements in GM: 
% 0.9875*10^5 ./ (1-exp(-2000/1240)) = 1.2333e+005
% The middle of these two is:
% Ingenia AMC = 2.1412*10^5

x.Q.BackGrSupprPulses       = 2; % '2' or '5' (this parameter will be ignored when M0='no_background_suppression').
% GE 3D FSE uses 5 pulses, Philips 2D & Siemens 3D use both 2 pulses
% This parameter is used to estimate decrease of labeling efficiency  (0.83 & 0.75 for 2 & 5 pulses respectively)
% When you have an M0, but no background suppression, this parameter should be set to '0'
x.readout_dim             = '2D'; % 2D or 3D
x.CLIP_ZERO               = 1; % zero clipping is better for registration but inconvenient if you have to correct control-label order. and create bias for estimates such as intra-scan variance (leave on)
x.SpikeRemovalThreshold 	= 0.05; % SNR needs to be improved by this rate before frames are excluded. 1 == disabling Spike Removal, 0.05==default
x.motion_correction       = 1; % 1 = on, 0 = off (leave on)
x.Quality                 = 1; % 1 = normal, 0 = low for fast try-out
x.DELETETEMP              = 0; % removes temporary files. saves lots of disk space, recommended.
x.DisableQuantification   = 0;% = 1 if the input images are already quantified CBF and you just want to run registration steps
x.Vendor                  = 'Philips';   % Options: 'GE_product' 'GE_WIP' 'Philips' 'Siemens', for applying vendor-specific scale factors
x.Sequence                = '2D_EPI' % options 3D_spiral, 3D_GRASE or 2D_EPI
x.Q.LabelingType            = 'CASL'; % Options: 'PASL' (pulsed Q2-TIPS) or 'CASL' (CASL/PCASL)
% NB: pulsed without Q2TIPS cannot be reliably quantified because the bolus width cannot be identified
% CASL & PCASL are both continuous ASL methods, identical quantification
x.Q.Initial_PLD            = 1800; % PLD, for 3D this is fixed for whole brain, for 2D this is the PLD of first acquired slice
x.Q.LabelingDuration              = 1800; % labeling duration
x.Q.SliceReadoutTime     = 'shortestTR';% for 2D, this is the time addition for each acquired slice later in time
%Do you have the protocol export (xml) from your PCASL? This gives you the minimal TR, from there you can calculate the SliceWise PLD addition.
% Other option = 'shortestTR'; % shortest TR enabled gives each sequence the minimal TR. This enables calculating slice delay per subject
% Other option = 'individual'; in this case the slice delay is specified in ASL4D_parms.mat (e.g. when using two different scanners)
% Examples:
% NOVICE: 37.06667
% AgeIV: 37.0625
% PreDiva baseline: 34.9;
% PreDiva crushed: 37.86667
% PreDiva non-crushed: 42.6 (probably other way around,concerning other studies, and the fact that crushing takes longer)
% GENFI Philips Bsup:   36.3;  % calculated by (minTR=3792-1525-1650)/17 slices
% GENFI Philips noBsup: 23.37; % although this probably differs slightly between sites
% Sleep study Oslo: = 39.91; 
% Cruise: = 42.1053 ms

%x.Q.LabelingEfficiency            = 0.85;            % labelling efficiency DEFAULT = 0.85
%x.Q.BloodT1                = 1650;            % T1 relaxation time of arterial blood DEFAULT=1650
x.Q.Lambda               = 0.9;              % Brain/blood water coefficient (mL 1H/ mL blood) or rho if in denominator of eq = 1/lambda = 1.1
%x.Q.T2art              = 165;              % T2* of arterial blood, used by Magdalena Sokolska UCL in simulations
x.Q.T2art               = 50;              % T2* of arterial blood, only used when no M0 image
x.Q.TissueT1                 = 1300; % 1300 measured with Turbo-QUASAR, default is 1240
%x.Q.TissueT1            = 1300; % T1 relaxation time of GM tissue DEFAULT=1240 ms @ 3T, 920 ms @ 1.5 T (Alsop MRM 2014)
%x.Q.ATT                 = 1500; % This is the arterial transit time

x.T1_DARTEL               = 0; % run T1  DARTEL
x.BILAT_FILTER            = 1; % 1 == run bilateral filter on raw EPI images, 2 == run filter on subtracted images
% This removes certain artifacts produced by fat suppression
% Bilateral filter 2 doesnt work that well for white matter because it uses
% WM to correct certain artefacts in the GM, so best try out the filter 1 

x.SkipIfNoFlair           = 1; % This will skip processing for any subject that doesn't contain a FLAIR
x.SkipIfNoASL             = 1; % This will skip processing for any subject that doesn't contain a ASL
x.SkipIfNoM0              = 1; % This will skip processing for any subject that doesn't contain a M0
% These settings determine whether we require certain files to be present, or otherwise to skip this subject
% By default we don't check this (because a population might not contain a
% FLAIR or M0, or this differs between subjects). This is useful when the
% data is incomplete, but one wants to start image processing nevertheless

%x.Q.nCompartments       = 2; % 1 = single-compartment quantification model, 2 = dual-compartment quantification model. 1 = default (concensus paper)

end
