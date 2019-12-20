function Ydiv = cat_vol_div(Ym,vx_vol,vx_volr)
% ----------------------------------------------------------------------
% Divergence helps to identify all gyri that should not be in the GM, but 
% helps to improve the WM. Divergence estimation is very memory intensive 
% so it is better to limit the resolution.
%
%   Ydiv = cat_vol_div(Ym[,vx_vol,vx_volr])
%  
%   Ym        .. input image
%   vx_vol    .. voxel resolution
%   vx_volr   .. lower voxel resolution for faster processing 
%                and smoother results
%
% ----------------------------------------------------------------------
% Robert Dahnke
% $Id: cat_vol_div.m 1310 2018-04-25 13:55:32Z dahnke $

  if nargin==0, help cat_vol_div; return; end
  if ~exist('vx_vol','var'), vx_vol = repmat(1.5,1,3); end % no reduction
  if numel(vx_vol)==1, vx_vol = repmat(vx_vol,1,3); end
  if ~exist('vx_volr','var'), vx_volr = min(1.5,vx_vol*3); end
  if numel(vx_volr)==1, vx_volr = repmat(vx_volr,1,3); end
  
  Ym = single(Ym); 
  Ynan = isnan(Ym); 
  
  % replace nan by neighbor values
  [D,I] = cat_vbdist(single(~Ynan),cat_vol_morph(~Ynan,'d',2)); Ym(D<2) = Ym(I(D<2)); 
  clear D I
  
  % low resolution ... % don't forget small animals...
  [Ymr,resT2] = cat_vol_resize(Ym,'reduceV',vx_vol,vx_volr,16,'cubic'); 
  clear Ym
  
  % gradients
  [gx,gy,gz]  = cat_vol_gradient3(max(1/3,Ymr)); clear Ymr
  
  % divergence function was too memory demanding for some systems
  [px,junk,junk] = cat_vol_gradient3(gx./resT2.vx_volr(1)); clear gx junk
  [junk,qy,junk] = cat_vol_gradient3(gy./resT2.vx_volr(2)); clear gy junk
  [junk,junk,rz] = cat_vol_gradient3(gz./resT2.vx_volr(3)); clear gz junk
  Ydivr = single(px) + single(qy) + single(rz); clear px qy rz
  
  % full resolution
  Ydiv  = cat_vol_resize(smooth3(Ydivr),'dereduceV',resT2,'cubic'); 
  Ydiv(Ynan) = 0;  % restore nan?
return