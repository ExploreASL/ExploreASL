function x = DATA_PAR( x )
%DATA_PAR Part of ExploreASL, loading basic study settings
% HJMM Mutsaerts 2018

% Define study
x.name               = 'Philips2DEPI';
x.subject_regexp     = '^(020|022|030|040)-\d{5}$';


% list parameters here
x.M0 				    = 'separate_scan';
x.Q.BackGrSupprPulses = 2;
x.readout_dim          = '2D'; % 2D or 3D
x.Quality              = 1; % 1 = normal, 0 = low for fast try-out

x.Vendor        	 	= 'Philips';   % Options: 'GE_product' 'GE_WIP' 'Philips' 'Siemens', for applying vendor-specific scale factors
x.Q.LabelingType      = 'CASL'; % Options: 'PASL' (pulsed Q2-TIPS) or 'CASL' (CASL/PCASL)
% NB: pulsed without Q2TIPS cannot be reliably quantified because the bolus width cannot be identified
% CASL & PCASL are both continuous ASL methods, identical quantification
x.Q.Initial_PLD         = 2025;
x.Q.LabelingDuration   = 1650;
x.Q.SliceReadoutTime   = 36.5278; % (minTR-PLD-labdur)/n_slices = (4990-1650-2025)/36 = 36.5278



end