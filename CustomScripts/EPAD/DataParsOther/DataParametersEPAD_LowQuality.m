function x = DataParameters( x )
%DataParameters EPAD ExploreASL
% HJMM Mutsaerts 2019

% Define study
x.name               = 'EPAD';
x.subject_regexp     = '^\d{3}EPAD\d*(|_\d*)$';
x.Quality  			 = false;
x.DELETETEMP 		 = true;

end