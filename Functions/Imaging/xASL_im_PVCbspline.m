function [imPVEC,imCBFrec,imResidual,FWHM] = xASL_im_PVCbspline(imCBF,imPV,bsplineNum)
%xASL_im_PVCbspline Partial volume correction (PVC) of ASL data using prior GM-,WM-partial volume maps using
% a global optimization with B-spline approximation of the correction weights rather than a local kernel.
%
% FORMAT:       [imPVEC,imCBFrec,imResidual,FWHM] = xASL_im_PVCbspline(imCBF,imPV[,bsplineNum])
%
% INPUT:
%   imCBF - uncorrected CBF image imCBF(X,Y,Z) (REQUIRED)
%   imPV  - PV-maps imPV(X,Y,Z,2) - WM/GM does not matter, you get the same
%           order in output as it is defined for input (REQUIRED)
%   bsplineNum - number of cubic B-splines for each dimension
%              - bsplineNum(I,J,K), it should be smaller than (X,Y,Z),
%                taking into consideration if the imCBF was upsampled
%              - (I,J,K) should stay the same for upsampled imCBF as for
%                the low resolution imCBF 
%                (OPTIONAL, DEFAULT = [15 15 1]
% OUTPUT: 
%   imPVEC - corrected CBF maps for both tissues - imPVEC(X,Y,Z,2)
%   imCBFrec - (X,Y,Z) - reconstruction of CBF from the corrected values
%                        and PV maps
%   imResidual(X,Y,Z) - difference between the reconstructed and original
%   FWHM - approximate Gaussian FWHM of the smoothness of the corrected CBF maps
%        - no that this is done in Pixels, so it needs to be multiplied by
%          the input image voxel size
% 
% ----------------------------------------------------------------------------------------------
% DESCRIPTION:  PVC of ASL data using prior GM-,WM-partial volume maps.
%               Follows the principles of the PVEc algorithm by I. Asllani (MRM, 2008).
%               The PV-corrected CBF_GM and CBF_WM maps are approximated using an
%               uniformly sampled cubic B-splines.
% ----------------------------------------------------------------------------------------------
% 
% REFERENCES:
% -Asllani I, Borogovac A, Brown TR. Regression algorithm correcting for 
%  partial volume effects in arterial spin labeling MRI. Magnetic Resonance 
%  in Medicine. 2008 Dec 1;60(6):1362-71.
% -Petr J, Mutsaerts HJ, De Vita E, Steketee RM, Smits M, Nederveen AJ, 
%  Hofheinz F, van den Hoff J, Asllani I. Effects of systematic partial 
%  volume errors on the estimation of gray matter cerebral blood flow with 
%  arterial spin labeling MRI. MAGMA 2018. DOI:10.1007/s10334-018-0691-y
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [imPVEC,imCBFrec,imResidual,FWHM] = xASL_im_PVCbspline(imCBF,imPV,[17 17 1])
% __________________________________
% Copyright 2015-2021 ExploreASL

if nargin<2
	error('CBF image and PV-maps have to be given on the input');
end

if nargin<3 || isempty(bsplineNum)
	bsplineNum = [15 15 1];
end

% Basic exclusion parameters
matchThresh = 0.0001;
bmatppThresh = 1;

% Checking the input
if sum(size(imCBF,1)<(3*bsplineNum(1)))
    error('B-spline number too high');
end

if sum(size(imCBF,2)<(3*bsplineNum(2)))
    error('B-spline number too high');
end

if (size(imCBF,3)==1) && (bsplineNum(3)>1)
    error('For single slice, bsplineNum(3) must be 1');
end

if (size(imCBF,3)<bsplineNum(3))
    error('The number of B-splines should not exceed the image size');
end

if (size(imCBF,3)>20) && (size(imCBF,3)<(1.2*bsplineNum(3)))
    error('B-spline number too high');
end

if max([size(imCBF),2] ~= size(imPV))
    error('The size of imCBF and imPV does not match');
end

% Sample the B-spline basis functions
splineX = xASL_im_cubicBspline(size(imCBF,1),bsplineNum(1));
splineY = xASL_im_cubicBspline(size(imCBF,2),bsplineNum(2));
splineZ = xASL_im_cubicBspline(size(imCBF,3),bsplineNum(3));

% Estimate the equivalent of Gaussian FWHM of the kernel
FWHM = 1.42*size(imCBF)./(bsplineNum-3);

Bmat = zeros([bsplineNum,2]);
BmatPP = zeros(size(Bmat));
fprintf('%s','B-spline calculation 1:   ');
for tt = 1:2
	for kk = 1:bsplineNum(3)
		resK = repmat(shiftdim(splineZ(:,kk)',-1),[size(imCBF,1),size(imCBF,2),1]);
		resK = resK.*imPV(:,:,:,tt);
		resKPV = sum(resK,3);
		resKCBF = sum(resK.*imCBF,3);
		for jj = 1:bsplineNum(2)
			resJ = repmat(splineY(:,jj)',[size(imCBF,1),1]);
			resJPV = sum(resKPV.*resJ,2);
			resJCBF = sum(resKCBF.*resJ,2);
			for ii = 1:bsplineNum(1)
				resI = splineX(:,ii);
				BmatPP(ii,jj,kk,tt) = sum(resJPV.*resI);
				Bmat(ii,jj,kk,tt) = sum(resJCBF.*resI);
			end
		end
		% keep track
		CurrentIt   = (tt-1)*bsplineNum(3) + kk-1;
		MaxIt       = bsplineNum(3)*2;
		xASL_TrackProgress(CurrentIt,MaxIt);
	end
end
Bmat = reshape(Bmat,[prod(bsplineNum)*2,1]);
fprintf('\n');
% BmatPP contains the sum of B-spline * pGM+pWM - this is later used to
% exclude the ill-conditioned variables from the inversion
BmatPP = reshape(BmatPP,[prod(bsplineNum)*2,1]);

% For these combinations of splineX*splineX, a zeros entry in Amat is
% produced, so they can be excluded from Amat calculation to speed things
% up
iMatch = zeros(bsplineNum(1),bsplineNum(1));
for i=1:bsplineNum(1)
	for ii=1:bsplineNum(1)
		iMatch(i,ii) = sum(splineX(:,i).*splineX(:,ii));
	end
end
jMatch = zeros(bsplineNum(2),bsplineNum(2));
for j=1:bsplineNum(2)
	for jj=1:bsplineNum(2)
		jMatch(j,jj) = sum(splineY(:,j).*splineY(:,jj));
	end
end
kMatch = zeros(bsplineNum(3),bsplineNum(3));
for k=1:bsplineNum(3)
    for kk=1:bsplineNum(3)
        kMatch(k,kk) = sum(splineZ(:,k).*splineZ(:,kk));
	end
end

Amat = zeros([bsplineNum,2,bsplineNum,2]);
fprintf('%s','B-spline calculation 2:   ');
for tt = 1:2
for t = tt:2
    resT = imPV(:,:,:,tt).*imPV(:,:,:,t);
    for kk = 1:bsplineNum(3)
    for k = kk:bsplineNum(3)
    if  (kMatch(kk,k)>matchThresh)
    	resK = (splineZ(:,kk)').*(splineZ(:,k)');
    	resK = repmat(shiftdim(resK,-1),[size(imCBF,1),size(imCBF,2),1]);
        resKT = sum(resK.*resT,3);
        for jj = 1:bsplineNum(2)
        for j = jj:bsplineNum(2)
        if (jMatch(jj,j)>matchThresh)   
            resJ = (splineY(:,jj)').*(splineY(:,j)');
            resJ = repmat(resJ,[size(imCBF,1),1]);
            resJKT = sum(resJ.*resKT,2);
            for ii = 1:bsplineNum(1)
                
            for i = ii:bsplineNum(1)
            if (iMatch(ii,i)>matchThresh)
                resI = splineX(:,ii).*splineX(:,i);
                val = sum(resJKT.*resI);

                Amat(ii,jj,kk,tt,i,j,k,t) = val;
                %Amat(ii,jj,kk,tt,i,j,k,t) = val;
                %Amat(ii,jj,kk,t,i,j,k,tt) = val;
                %Amat(ii,jj,k,tt,i,j,kk,t) = val;
                %Amat(ii,jj,k,t,i,j,kk,tt) = val;
                %Amat(ii,j,kk,tt,i,jj,k,t) = val;
                %Amat(ii,j,kk,t,i,jj,k,tt) = val;
                %Amat(ii,j,k,tt,i,jj,kk,t) = val;
                %Amat(ii,j,k,t,i,jj,kk,tt) = val;
                %Amat(i,jj,kk,tt,ii,j,k,t) = val;
                %Amat(i,jj,kk,t,ii,j,k,tt) = val;
                %Amat(i,jj,k,tt,ii,j,kk,t) = val;
                %Amat(i,jj,k,t,ii,j,kk,tt) = val;
                %Amat(i,j,kk,tt,ii,jj,k,t) = val;
                %Amat(i,j,kk,t,ii,jj,k,tt) = val;
                %Amat(i,j,k,tt,ii,jj,kk,t) = val;
                %Amat(i,j,k,t,ii,jj,kk,tt) = val;
			end
			end
			end
		end
		end
		end
	end
	end
	% keep track
	CurrentIt   = (tt-1)*2*bsplineNum(3) + (t-1)*bsplineNum(3) + kk;
	MaxIt       = bsplineNum(3)*4;
	xASL_TrackProgress(CurrentIt,MaxIt);
	end
end
end
fprintf('\n');
% Ensure the full symmetricity of Amat in I,J,K,T
for tt = 1:2
    for t = tt:2
        Amat(:,:,:,t,:,:,:,tt) = Amat(:,:,:,tt,:,:,:,t);
	end
end

for kk = 1:bsplineNum(3)
    for k = kk:bsplineNum(3)
        Amat(:,:,k,:,:,:,kk,:) = Amat(:,:,kk,:,:,:,k,:);
	end
end
for jj = 1:bsplineNum(2)
    for j = jj:bsplineNum(2)
        Amat(:,j,:,:,:,jj,:,:) = Amat(:,jj,:,:,:,j,:,:);
	end
end
for ii = 1:bsplineNum(1)
    for i = ii:bsplineNum(1)
        Amat(i,:,:,:,ii,:,:,:) = Amat(ii,:,:,:,i,:,:,:);
	end
end
Amat = reshape(Amat,[prod(bsplineNum)*2,prod(bsplineNum)*2]);           

% For these B-splines, no solution exists so they are set to zero and
% excluded from the solution
indPP = find(BmatPP>bmatppThresh);

% Matrix inversion
AmatInv = inv(Amat(indPP,indPP));
sol = AmatInv*Bmat(indPP);

% Solution vector
G = zeros(size(Bmat));
% Full solution vector
G(indPP) = sol;

G = reshape(G,[bsplineNum,2]);
       
% Recalculate the CBF-GM maps from the B-spline parameters
resX = splineX*reshape(G,[bsplineNum(1),bsplineNum(2)*bsplineNum(3)*2]);
resX = reshape(resX,[size(imCBF,1),bsplineNum(2),bsplineNum(3),2]);
resY = zeros([size(imCBF,1),size(imCBF,2),bsplineNum(3),2]);
for x=1:size(imCBF,1)
    res = splineY*reshape(resX(x,:,:,:),[bsplineNum(2),bsplineNum(3)*2]);
    resY(x,:,:,:) = reshape(res,[1,size(imCBF,2),bsplineNum(3),2]);
end

resZ = zeros([size(imCBF),2]);
for x=1:size(imCBF,1)
    for y=1:size(imCBF,2)
        res = splineZ*squeeze(resY(x,y,:,:));
        resZ(x,y,:,:) = reshape(res,[1,1,size(imCBF,3),2]);
	end
end

% Calculate the extra outputs
imPVEC = resZ;
imCBFrec = imPVEC(:,:,:,1).*imPV(:,:,:,1) + imPVEC(:,:,:,2).*imPV(:,:,:,2);
imResidual = abs(imCBF-imCBFrec);

% TODO:
% NULL the negative values the negative values somehow
% variable output
% GNU license
% It is possitive definite and symmetric - use that for Amat
% Use Cholesky decomposition solution
% Remove them (those not necessary) already for Amat calculation
% Avoid negative numbers in the solution
% Commentaries + clean the code
% Big numbers outside of brain
return;
% 
% The cubic B-spline calculation
function spline = xASL_im_cubicBspline(N,bsplineNum)

H = N./(bsplineNum-3);

spline = repmat(1:bsplineNum,[N,1]);
X = repmat((0.5:(N-0.5))',[1,bsplineNum]);
spline = X/H - spline + 2;

ind0 = find(abs(spline)<1);
ind1 = find((abs(spline)<2).*(abs(spline)>=1));
ind2 = abs(spline)>=2;

spline(ind2) = 0;
spline(ind1) = ((2-abs(spline(ind1))).^3)/6;
spline(ind0) = 2/3 - (spline(ind0).^2) + abs(spline(ind0).^3)/2;

return;
