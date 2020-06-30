function P = cat_main_clean_gwc(P,level,new)
%==========================================================================
% Cleanup function that removes non-brain tissue (e.g. meninges) by means
% of morphological operations to cleanup WM, GM, and CSF tissue.
% Successor of the cg_cleanup_gwc function of VBM8. 
% Include a new brain limitation that remove/add empty space around the
% brain for speedup (no functional difference). 
% Moreover, a new morpholocial cleanup close to the skull was added to
% remove larger unwanted parts of head tissue that is used for the kamap
% preprocessing pipeline (see cat_main_kamap, 201812).
%
%  function P = cat_main_clean_gwc(P[,level,new])
% 
%  P     .. 4D uint8 matrix of tissue classes GM, WM, and CSF
%  level .. controls strength of corrections (values=[1,2]; default=1)
%  new   .. use new additional cleanup (default=0)
%__________________________________________________________________________
% Christian Gaser, Robert Dahnke
% $Id: cat_main_clean_gwc.m 1401 2018-12-04 22:36:50Z dahnke $

if nargin<2, level  = 1; end
if nargin<3, new    = 0; end

% remove empty space for speed up
sP = sum(P,4); 
for i=1:size(P,4), [P2(:,:,:,i),BB] = cat_vol_resize(P(:,:,:,i),'reduceBrain',ones(1,3),2,sP); end  %#ok<AGROW>
P = P2; clear sP P2;

% New additional harder cleanup close to the skull to remove meninges.
% Added due to problems with the alternative cat_main_kamap segmentation. 
% However, this should also help in other cases and should not create to 
% large problems. >> TEST IT! (RD: 201812)
%--------------------------------------------------------------------------
if new
  Yp0  = single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255*2 + single(P(:,:,:,2))/255*3;
  Ybe  = cat_vol_morph(cat_vol_morph( cat_vol_morph(Yp0>0.5,'ldc',1), 'de', 5 ), 'lc'); % area close to the skull
  % mask WM
  Ymsk = Ybe | cat_vol_morph( cat_vol_morph( Yp0>2.5 , 'do' , 0.5 + level/2 ) , 'l' , [10 0.1]); 
  P(:,:,:,2) = cat_vol_ctype(single(P(:,:,:,2)) .* cat_vol_smooth3X(cat_vol_morph(Ymsk,'dd',1.5),0.5)); 
  Yp0  = single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255*2 + single(P(:,:,:,2))/255*3; % update
  % mask GM (iter
  Ywd  = cat_vol_morph( Yp0>2.5 , 'dd' , 5 - level); % area close to the WM
  for i=1:2
    Ymsk = cat_vol_morph(Ybe | (Ywd & cat_vol_morph( Yp0>1.5 , 'do' , level/2+1.5 )),'ldc',1); 
    P(:,:,:,1) = cat_vol_ctype(single(P(:,:,:,1)) .* cat_vol_smooth3X(cat_vol_morph(Ymsk,'dd',1.5),0.5)); 
  end
  Yp0  = single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255*2 + single(P(:,:,:,2))/255*3; % update
  % CSF
  for i=1:2
    Ymsk = cat_vol_morph(Ybe | (cat_vol_morph( Yp0>0.95 , 'do' , level+2.5 )),'ldc',1); 
    P(:,:,:,3) = cat_vol_ctype(single(P(:,:,:,3)) .* cat_vol_smooth3X(cat_vol_morph(Ymsk,'dd',1.5),0.5)); 
  end
  clear Yp0 Ybe Ymsk Ywd;
end

%b = P(:,:,:,2);
b = cat_vol_morph(P(:,:,:,2)>128,'l',[10 0.1])*255;

% Build a 3x3x3 seperable smoothing kernel
%--------------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

th1 = 0.15;
if level>1, th1 = 0.2; end
% Erosions and conditional dilations
%--------------------------------------------------------------------------
niter  = 32;
niter2 = 32;
spm_progress_bar('Init',niter+niter2,'Extracting Brain','Iterations completed');
fprintf('Clean Up 1:   ');
for j=1:niter
    xASL_TrackProgress(j,niter);
    if j>2, th=th1; else th=0.6; end  % Dilate after two its of erosion
    for i=1:size(b,3)
        gp       = single(P(:,:,i,1));
        wp       = single(P(:,:,i,2));
        bp       = single(b(:,:,i))/255;
        bp       = (bp>th).*(wp+gp);
        b(:,:,i) = cat_vol_ctype(round(bp));
    end
    spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
    spm_progress_bar('Set',j);
end

% Also clean up the CSF.
fprintf('\n');
fprintf('Also clean up the CSF:   ');
if niter2 > 0,
    c = b;
    for j=1:niter2
        xASL_TrackProgress(j,niter2);
        for i=1:size(b,3)
            gp       = single(P(:,:,i,1));
            wp       = single(P(:,:,i,2));
            cp       = single(P(:,:,i,3));
            bp       = single(c(:,:,i))/255;
            bp       = (bp>th).*(wp+gp+cp);
            c(:,:,i) = cat_vol_ctype(round(bp));
        end
        spm_conv_vol(c,c,kx,ky,kz,-[1 1 1]);
        spm_progress_bar('Set',j+niter);
    end
end

th = 0.05;
fprintf('\n');
fprintf('Some slice-work:   ');
for i=1:size(b,3)
    xASL_TrackProgress(i,size(b,3));
    slices = cell(1,size(P,4));
    for k1=1:size(P,4),
        slices{k1} = single(P(:,:,i,k1))/255;
    end
    bp        = single(b(:,:,i))/255;
    bp        = ((bp>th).*(slices{1}+slices{2}))>th;
    slices{1} = slices{1}.*bp;
    slices{2} = slices{2}.*bp;

    if niter2>0,
        cp        = single(c(:,:,i))/255;
        cp        = ((cp>th).*(slices{1}+slices{2}+slices{3}))>th;
        slices{3} = slices{3}.*cp;
    end
    if numel(slices)>=5
      slices{5} = slices{5}+1e-4; % Add something to the soft tissue class
    end
    tot       = zeros(size(bp))+eps;
    for k1=1:size(P,4),
        tot   = tot + slices{k1};
    end
    for k1=1:size(P,4),
        P(:,:,i,k1) = cat_vol_ctype(round(slices{k1}./tot*255));
    end 
end

% add previously removed empty space 
for i=1:size(P,4), P2(:,:,:,i) = cat_vol_resize(P(:,:,:,i),'dereduceBrain',BB); end; 
P = P2; clear P2; 

spm_progress_bar('Clear');
return;
%==========================================================================
