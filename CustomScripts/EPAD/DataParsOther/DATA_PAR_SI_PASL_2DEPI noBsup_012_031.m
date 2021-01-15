function x = DATA_PAR( x )
%DATA_PAR Part of ExploreASL, loading basic study settings
% HJMM Mutsaerts 2018

% Define study
x.name               = 'Siemens2DEPI';
x.subject_regexp     = '^(012|031)-\d{5}$';


% list parameters here
x.M0 				    = 'UseControlAsM0';  % Dennis obtained in PET-MRI study, rescale Ingenia-Intera later
x.Q.BackgroundSuppressionNumberPulses = 0;
x.readout_dim          = '2D'; % 2D or 3D
x.Quality              = 1; % 1 = normal, 0 = low for fast try-out

x.Vendor        	 	= 'Siemens';   % Options: 'GE_product' 'GE_WIP' 'Philips' 'Siemens', for applying vendor-specific scale factors
x.Q.LabelingType      = 'PASL'; % Options: 'PASL' (pulsed Q2-TIPS) or 'CASL' (CASL/PCASL)
% NB: pulsed without Q2TIPS cannot be reliably quantified because the bolus width cannot be identified
% CASL & PCASL are both continuous ASL methods, identical quantification
x.Q.Initial_PLD         = 2000; % Philips Utrecht == instellingen origineel van Thijs (v??r white paper)
x.Q.LabelingDuration   = 800; %?
x.Q.SliceReadoutTime   = 35; % ? (minTR-PLD-labdur)/n_slices = (5000-2000)/36 = 83.33



end