function x = DATA_PAR(x)
%DATA_PAR Part of ASL pipeline, loading basic study settings
% AMC ASL pipeline, HJMM Mutsaerts 2014


% Define study
x.name               = 'RADAR_Bristol';
x.ROOT               = cd;
x.subject_regexp     = '^\d{7}_\d$';


% list parameters here
x.M0 				    = 10^4;  % Dennis obtained in PET-MRI study, rescale Ingenia-Intera later
x.Q.BackGrSupprPulses 	= 2;
x.readout_dim          = '3D'; % 2D or 3D

x.Vendor        	 	= 'Siemens';   % Options: 'GE_product' 'GE_WIP' 'Philips' 'Siemens', for applying vendor-specific scale factors
x.Q.LabelingType        	= 'PASL'; % Options: 'PASL' (pulsed Q2-TIPS) or 'CASL' (CASL/PCASL)
% NB: pulsed without Q2TIPS cannot be reliably quantified because the bolus width cannot be identified
% CASL & PCASL are both continuous ASL methods, identical quantification
x.Q.Initial_PLD       = 2000;  %  % Philips Utrecht == instellingen origineel van Thijs (v??r white paper)
x.Q.LabelingDuration  = 800; % (90rfblocks x 18.4ms). previously processed with 1650, hardly any difference
x.Q.SliceReadoutTime  = 0; % minTR = (3980 - 1656 - 1525)/16  = 49.9375. Previous educated guess was 51.5625, hardly different


end