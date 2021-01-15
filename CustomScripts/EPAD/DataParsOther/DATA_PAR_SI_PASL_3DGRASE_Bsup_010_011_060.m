function x = DATA_PAR( x )
%DATA_PAR Part of ExploreASL, loading basic study settings
% HJMM Mutsaerts 2018

% Define study
x.name               = 'Siemens3DGRASE';
x.subject_regexp     = '^(010|011|060)-\d{5}$';



% list parameters here
x.M0 				    = 'UseControlAsM0'; % this has background suppression, but doesn't have an M0, so let's divide by control image
x.Q.BackgroundSuppressionNumberPulses = 2;
x.readout_dim          = '3D'; % 2D or 3D
x.Quality              = 1; % 1 = normal, 0 = low for fast try-out

x.Vendor        	 	= 'Siemens';   % Options: 'GE_product' 'GE_WIP' 'Philips' 'Siemens', for applying vendor-specific scale factors
x.Q.LabelingType      = 'PASL'; % Options: 'PASL' (pulsed Q2-TIPS) or 'CASL' (CASL/PCASL)
% NB: pulsed without Q2TIPS cannot be reliably quantified because the bolus width cannot be identified
% CASL & PCASL are both continuous ASL methods, identical quantification
x.Q.Initial_PLD         = 2000;
x.Q.LabelingDuration   = 800;
x.Q.SliceReadoutTime   = 0;
x.Q.NumberSegments 	   = 4;
x.Q.Sequence 		   = '3D_GRASE';
x.Q.NumberOfAverages 	= 10;

end