function [x] = xASL_adm_DefineASLSequence(x)
%xASL_adm_DefineASLSequence Obtain ASL sequence type for sequence-specific
%image processing
%
% FORMAT: [x] = xASL_adm_DefineASLSequence(x)
%
% INPUT:
%   x               - x structure containing all input parameters (REQUIRED)
%   x.readout_dim   - dimensionality of readout (2D or 3D)
%   x.Vendor        - Either 'GE', 'Philips', 'Siemens'
% OUTPUT:
%   x               - x structure containing all output parameters
%   x.Sequence      - sequence type (readout)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This ExploreASL function tries to check what ASL sequence is
% being processed, if this was not already defined in x.Sequence.
% It does so by checking known combinations of readout dimensionality
% (x.readout_dim) and vendor, knowing the product sequences of the vendors.
%
% EXAMPLE: x = xASL_adm_DefineASLSequence(x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL




% Obtain ASL sequence
if ~isfield(x,'Sequence') && isfield(x,'readout_dim')
    if strcmpi(x.readout_dim,'2D')
       x.Sequence = '2D_EPI'; % assume that 2D is 2D EPI, irrespective of vendor
    elseif strcmpi(x.readout_dim,'3D') && ( ~isempty(regexpi(x.Vendor,'Philips')) || ~isempty(regexpi(x.Vendor,'Siemens')) )
           x.Sequence = '3D_GRASE'; % assume that 3D Philips or Siemens is 3D GRASE
    elseif strcmpi(x.readout_dim,'3D') && ~isempty(regexpi(x.Vendor,'GE'))
           x.Sequence = '3D_spiral'; % assume that 3D GE is 3D spiral
    end
elseif ~isfield(x,'Sequence') && ~isfield(x,'readout_dim')
    warning('No x.Sequence defined');
    fprintf('If there are multiple sequence types, this needs to be implemented yet here\n');
    fprintf('Otherwise, please define x.Sequence\n');
    fprintf('Setting x.Sequence=3D spiral to fool the pipeline here\n');
    fprintf('As this sequence doesnt have any susceptibility masks\n');
    fprintf('Note that this disables any masking of susceptibility signal dropout areas\n');
    x.Sequence = '3D_spiral';
    x.readout_dim = '3D';
end




end
