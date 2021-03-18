function xASL_wrp_DARTELSaveIntermedTrans(Yy, u, odim, rdim, idim, Mar, mat, M0, M1, nameOut, numIteration)
%xASL_wrp_DARTELSaveIntermedTrans Saves the intermediate results, the transformation field of DARTEL
%
% FORMAT: xASL_wrp_DARTELSaveIntermedTrans(Yy, u, odim, rdim, idim, Mar, mat, M0, M1, nameOut, numIteration)
%
% INPUT:
% Yy           - internal parameter of the CAT12 - segmentation (REQUIRED)
% u            - internal parameter of the CAT12 - segmentation (REQUIRED)
% odim         - internal parameter of the CAT12 - segmentation (REQUIRED)
% rdim         - internal parameter of the CAT12 - segmentation (REQUIRED)
% idim         - internal parameter of the CAT12 - segmentation (REQUIRED)
% Mar          - internal parameter of the CAT12 - segmentation (REQUIRED)
% mat          - internal parameter of the CAT12 - segmentation (REQUIRED)
% M0           - internal parameter of the CAT12 - segmentation (REQUIRED)
% M1           - internal parameter of the CAT12 - segmentation (REQUIRED)
% nameOut      - name of the output file (REQUIRED)
% numIteration - number of the registration iteration (REQUIRED)
%
% OUTPUT: Saves the deformation field under a filename ./mri/y_nameOut_numIteration.nii
%
% DESCRIPTION: This function is called from the CAT12 segmentation function to save the intermediate results
%              of the DARTEL transformation. Normally, the registration only saves the final results - 
%              the final transformation field. This function enables to save also the intermediate transformation field.
%              It takes all the internal variables from the transformation and save the field to a sub-directory 'mri'
%              that normally contains all the intermediate results of the CAT12 segmentation. It adds a 'y_' prefix
%              and adds the specified iteration number as postfix.
%
% EXAMPLE: xASL_wrp_DARTELSaveIntermedTrans(Yy, u, odim, rdim, idim, Mar, mat, M0, M1, 'T1',3)
%          Saves a file in ./mri/y_T1_3.nii
% __________________________________
% Copyright 2015-2021 ExploreASL

if nargin <11
	error('Requires 11 input parameters');
end

%% Calculates the transformation field
y0      = spm_dartel_integrate(reshape(u,[rdim(1:3) 1 3]),[0 1], 6);
prm     = [3 3 3 0 0 0];
Coef    = cell(1,3);
Coef{1} = spm_bsplinc(y0(:,:,:,1),prm);
Coef{2} = spm_bsplinc(y0(:,:,:,2),prm);
Coef{3} = spm_bsplinc(y0(:,:,:,3),prm);
clear y0;

[t1,t2] = ndgrid(1:idim(1),1:idim(2),1); 
t3     = 1:idim(3);

Yyd = Yy;
for z=1:idim(3)
	[t11,t22,t33] = defs2(Coef,z,Mar,prm,t1,t2,t3);
	Yyd(:,:,z,1) = t11;
	Yyd(:,:,z,2) = t22;
	Yyd(:,:,z,3) = t33;
end
clear Coef t1 t2 t3 t11 t22 t33 z

M = mat\M1;
for i=1:size(Yyd,3)
	t1          = Yyd(:,:,i,1);
	t2          = Yyd(:,:,i,2);
	t3          = Yyd(:,:,i,3);
	Yyd(:,:,i,1) = M(1,1)*t1 + M(1,2)*t2 + M(1,3)*t3 + M(1,4);
	Yyd(:,:,i,2) = M(2,1)*t1 + M(2,2)*t2 + M(2,3)*t3 + M(2,4);
	Yyd(:,:,i,3) = M(3,1)*t1 + M(3,2)*t2 + M(3,3)*t3 + M(3,4);
end
clear t1 t2 t3 M;

vx_vols  = sqrt(sum(M0(1:3,1:3).^2));
vx_volt  = sqrt(sum(M1(1:3,1:3).^2));
interpol = any(vx_vols>vx_volt) + any(vx_vols/2>vx_volt);
if interpol
	Yx = zeros(min([inf inf inf 3],size(Yyd)*2^interpol - 2^interpol + 1),'single');
	if interpol
		for i=1:3, Yx(:,:,:,i) = single(interp3(Yyd(:,:,:,i),interpol,'cubic')); end % linear
	end
else
	Yx = Yyd;
end

clear vx_vols vx_volt Yyd interpol;	      

%% Saves the file
% Creates the output pathname
[Fpath,Fname] = spm_fileparts(nameOut);
nameOut  = fullfile(Fpath,'mri',['y_' Fname '_' num2str(numIteration) '.nii']);
fprintf('%s %d %s %s\n','Saving DARTEL intermediate step',numIteration,'to a file:',nameOut);

% Creates and saves the NIFTI file
Yy2       = spm_diffeo('invdef',Yx,odim,eye(4),M0);
N         = nifti;
N.dat     = file_array(nameOut,[odim(1:3),1,3],'float32',0,1,0);
N.mat     = M1;
N.mat0    = M1;
N.descrip = 'Deformation';
create(N);
N.dat(:,:,:,:,:) = reshape(Yy2,[odim,1,3]);
end

%=======================================================================
function [x1,y1,z1] = defs2(sol,z,M,prm,x0,y0,z0)
  iM = inv(M);
  z01 = z0(z)*ones(size(x0));

  x1a  = iM(1,1)*x0 + iM(1,2)*y0 + iM(1,3)*z01 + iM(1,4);
  y1a  = iM(2,1)*x0 + iM(2,2)*y0 + iM(2,3)*z01 + iM(2,4);
  z1a  = iM(3,1)*x0 + iM(3,2)*y0 + iM(3,3)*z01 + iM(3,4);

  x1 = spm_bsplins(sol{1},x1a,y1a,z1a,prm);
  y1 = spm_bsplins(sol{2},x1a,y1a,z1a,prm);
  z1 = spm_bsplins(sol{3},x1a,y1a,z1a,prm);
end         
          

    
