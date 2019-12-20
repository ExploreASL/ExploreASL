function x = DATA_PAR( x )
%DATA_PAR Part of ExploreASL, loading basic study settings
% HJMM Mutsaerts 2018

% Define study
x.name               = 'GE_3Dspiral_VUmcData';
x.subject_regexp     = '^060-\d{5}$';


% list parameters here
x.M0 				    = 'separate_scan';
x.Q.BackGrSupprPulses = 5;
x.readout_dim          = '3D'; % 2D or 3D
x.Quality              = 1; % 1 = normal, 0 = low for fast try-out

x.Vendor        	 	= 'GE_WIP';   % Options: 'GE_product' 'GE_WIP' 'Philips' 'Siemens', for applying vendor-specific scale factors
x.Q.LabelingType      = 'CASL'; % Options: 'PASL' (pulsed Q2-TIPS) or 'CASL' (CASL/PCASL)
% NB: pulsed without Q2TIPS cannot be reliably quantified because the bolus width cannot be identified
% CASL & PCASL are both continuous ASL methods, identical quantification
x.Q.Initial_PLD         = 1525; % Philips Utrecht == instellingen origineel van Thijs (v??r white paper)
x.Q.LabelingDuration   = 1450; %?

end