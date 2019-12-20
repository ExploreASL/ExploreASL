function x = DATA_PAR( x )
%DATA_PAR Part of ASL pipeline, loading basic study settings
% AMC ASL pipeline, HJMM Mutsaerts 2014

% list parameters here
x.subject_regexp     = '^MCI-\d{4}$';

x.M0 				    = 'separate_scan';  % Dennis obtained in PET-MRI study, rescale Ingenia-Intera later
x.Q.BackGrSupprPulses 	= 4;
x.readout_dim          = '3D'; % 2D or 3D
x.QUALITY              = 1; % 1 = normal, 0 = low for fast try-out
x.DELETETEMP           = 1; % removes temporary files. saves lots of disk space, recommended.

x.Vendor        	 	= 'Philips';   % Options: 'GE_product' 'GE_WIP' 'Philips' 'Siemens', for applying vendor-specific scale factors
x.Q.LabelingType        	= 'CASL'; % Options: 'PASL' (pulsed Q2-TIPS) or 'CASL' (CASL/PCASL)
% NB: pulsed without Q2TIPS cannot be reliably quantified because the bolus width cannot be identified
% CASL & PCASL are both continuous ASL methods, identical quantification
x.Q.Initial_PLD         = 1800; % Philips Utrecht == instellingen origineel van Thijs (vóór white paper)
x.Q.LabelingDuration           = 1800; %?
x.Q.SliceReadoutTime  = 0; % (minTR-PLD-labdur)/n_slices = (3919-1650-1525)/17 = 

end