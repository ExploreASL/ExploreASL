function Ycls = cat_main_updateWMHs1585(Ym,Ycls,Yy,tpm,job,res,trans)
% ______________________________________________________________________
% Handle WMHs in the segmentation.
%  
%   Ycls = cat_main_updateWMHs1585(Ym,Ycls,Yy,job,trans)
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_main_updateWMHs.m 1393 2018-11-19 08:36:43Z dahnke $

  % update WMHs?
  if numel(Ycls)>6
    if isfield(trans,'warped')
      %% load template 
      VwmA = spm_vol([job.extopts.templates{end},',2']);  
      
      if any( VwmA.dim ~= trans.warped.odim )
        % interpolation
        yn = numel(trans.warped.y); 
        p  = ones([4,yn/3],'single'); 
        p(1,:) = trans.warped.y(1:yn/3);
        p(2,:) = trans.warped.y(yn/3+1:yn/3*2);
        p(3,:) = trans.warped.y(yn/3*2+1:yn);
        amat   = VwmA.mat \ trans.warped.M1; 
        p      = amat(1:3,:) * p;

        Yy = zeros([res.image(1).dim(1:3),3],'single'); 
        Yy(1:yn/3)        = p(1,:);
        Yy(yn/3+1:yn/3*2) = p(2,:);
        Yy(yn/3*2+1:yn)   = p(3,:);

        Yy = double(Yy); 
      else
        Yy = double(trans.warped.y);
      end
      YwmA = single( spm_sample_vol( VwmA ,Yy(:,:,:,1),Yy(:,:,:,2),Yy(:,:,:,3),1)); YwmA = reshape(YwmA,size(Ym)); 
    else
      %% load template 
      VwmA = spm_vol([job.extopts.templates{end},',2']);  
      
      if any( VwmA.dim ~= size(Yy(:,:,:,1)) )
        % interpolation
        yn = numel(trans.atlas.Yy); 
        p  = ones([4,yn/3],'single'); 
        p(1,:) = trans.atlas.Yy(1:yn/3);
        p(2,:) = trans.atlas.Yy(yn/3+1:yn/3*2);
        p(3,:) = trans.atlas.Yy(yn/3*2+1:yn);
        amat   = VwmA.mat \ tpm.M; 
        p      = amat(1:3,:) * p;

        Yy = zeros([res.image(1).dim(1:3),3],'single'); 
        Yy(1:yn/3)        = p(1,:);
        Yy(yn/3+1:yn/3*2) = p(2,:);
        Yy(yn/3*2+1:yn)   = p(3,:);

        Yy = double(Yy); 
      else
        Yy = double(trans.warped.y);
      end
      YwmA = single( spm_sample_vol( VwmA ,Yy(:,:,:,1),Yy(:,:,:,2),Yy(:,:,:,3),1)); YwmA = reshape(YwmA,size(Ym)); 
    end

    % transfer tissue from GM to WMH class
    Yclst   = cat_vol_ctype( single(Ycls{7} ) .* YwmA ); 
    Ycls{2} = Ycls{2} + (Ycls{7} - Yclst);
    Ycls{7} = Yclst;

  end

  
  if job.extopts.WMHC<2 && numel(Ycls)>6
    Ycls{1} = Ycls{1} + Ycls{7}; % WMH as GM
  elseif job.extopts.WMHC==2
    Ycls{2} = Ycls{2} + Ycls{7}; % WMH as WM 
  elseif job.extopts.WMHC>=3
    Ycls{2} = Ycls{2} + Ycls{7}; % WMH as own class
  end

end