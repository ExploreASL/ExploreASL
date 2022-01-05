function [IM] = xASL_im_ZeroEdges(IM, EdgeThicknessPerc)
%xASL_im_ZeroEdges Set the edges of the image to zero
%
% FORMAT: [IM] = xASL_im_ZeroEdges(IM[, EdgeThicknessPerc])
%
% INPUT:
%   IM - image matrix to zero edges for
%   EdgeThicknessPerc  - thickness of edges relative to image matrix size (%) (OPTIONAL, DEFAULT = 0.05)
%
% OUTPUT: IM - image matrix for which edges have been set to zero
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Resampling can sometimes give some strange errors near image edges. These should be NaNs,
% but sometimes can be zeros or ones, or even weird numbers. For resampling, NaNs should be set to 0 (this is done
% in another function) as they can influence the resampling (depending on the transformation matrix). To be sure
% that the edges are nicely fixed, this function sets a border at the image matrix edges to zero.
%
% EXAMPLE: ImOut = xASL_im_ZeroEdges(ImIn);
% __________________________________
% Copyright 2015-2019 ExploreASL
%
% 2019-05-02 HJM


if nargin<2 || isempty(EdgeThicknessPerc)
    EdgeThicknessPerc = 0.05;
end

% Get thickness edge
nVox = round((size(IM)./2).*EdgeThicknessPerc);
% Correct edges (set to 0)
IM(1:nVox(1)      ,:              ,:              ) = 0;
IM(end-nVox(1):end,:              ,:              ) = 0;
IM(:              ,1:nVox(2)      ,:              ) = 0;
IM(:              ,end-nVox(2):end,:              ) = 0;
IM(:              ,:              ,1:nVox(3)      ) = 0;
IM(:              ,:              ,end-nVox(3):end) = 0;


end

