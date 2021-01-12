% +BIDS - The bids-matlab library
%
% Contents
%   layout   - Parse a directory structure formated according to the BIDS standard
%   query    - Query a directory structure formated according to the BIDS standard
%   report   - Create a short summary of the acquisition parameters of a BIDS dataset
%   validate - BIDS Validator
%
%   util.jsondecode - Decode JSON-formatted file
%   util.jsonencode - Encode JSON-formatted file
%   util.tsvread    - Load text and numeric data from tab-separated-value file
%   util.tsvwrite   - Save text and numeric data to tab-separated-value file
% __________________________________________________________________________
%
% BIDS (Brain Imaging Data Structure): https://bids.neuroimaging.io/
%   The brain imaging data structure, a format for organizing and
%   describing outputs of neuroimaging experiments.
%   K. J. Gorgolewski et al, Scientific Data, 2016.
% __________________________________________________________________________
%
% BIDS-MATLAB is a library that aims at centralising MATLAB/Octave tools
% for interacting with datasets conforming to the BIDS format.
% See https://github.com/bids-standard/bids-matlab for more details.
