function xASL_wrp_GSSaveIntermedTrans(y, idim, odim, rdim, M0, M1, R, M1t, M1r, nameOut, numIteration)
%xASL_wrp_GSSaveIntermedTrans Saves the intermediate results, the transformation field of the Geodesic shooting
%
% FORMAT: xASL_wrp_GSSaveIntermedTrans(y, idim, odim, rdim, M0, M1, R, M1t, M1r, nameOut, numIteration)
%
% INPUT:
% y            - internal parameter of the CAT12 - segmentation (REQUIRED)
% idim         - internal parameter of the CAT12 - segmentation (REQUIRED)
% odim         - internal parameter of the CAT12 - segmentation (REQUIRED)
% rdim         - internal parameter of the CAT12 - segmentation (REQUIRED)
% M0           - internal parameter of the CAT12 - segmentation (REQUIRED)
% M1           - internal parameter of the CAT12 - segmentation (REQUIRED)
% R            - internal parameter of the CAT12 - segmentation (REQUIRED)
% M1t          - internal parameter of the CAT12 - segmentation (REQUIRED)
% M1r          - internal parameter of the CAT12 - segmentation (REQUIRED)
% nameOut      - name of the output file (REQUIRED)
% numIteration - number of the registration iteration (REQUIRED)
%
% OUTPUT: Saves the deformation field under a filename ./mri/y_nameOut_numIteration.nii
%
% DESCRIPTION: This function is called from the CAT12 segmentation function to save the intermediate results
%              of the Geodesic shooting transformation. Normally, the registration only saves the final results - 
%              the final transformation field. This function enables to save also the intermediate transformation field.
%              It takes all the internal variables from the transformation and save the field to a sub-directory 'mri'
%              that normally contains all the intermediate results of the CAT12 segmentation. It adds a 'y_' prefix
%              and adds the specified iteration number as postfix.
%
% EXAMPLE: xASL_wrp_GSSaveIntermedTrans(y,idim,odim,rdim,M0,M1,R,M1t,M1r,'T1',3)
%          Saves a file in ./mri/y_T1_3.nii
% __________________________________
% Copyright 2015-2021 ExploreASL

if nargin <11
	error('Requires 11 input parameters');
end

%% Calculates the transformation field
% Update for output resolution
if any(odim ~= rdim)
	eyev = eye(4); eyev(1:end-1) = eyev(1:end-1) * M1t(1)./M1r(1);
	yid  = zeros([odim 3],'single');
	for k1=1:3
		for i=1:odim(3)
			yid(:,:,i,k1) = single(spm_slice_vol(y(:,:,:,k1),eyev*spm_matrix([0 0 i]),odim(1:2),[1,NaN])) / eyev(1); % adapt for res
		end
	end
else
	yid = y;
end
         
yi  = spm_diffeo('invdef',yid,idim,inv(M1t\R*M0),eye(4));           % output yi in anatomical resolution
clear eyev yid;

% interpolation for improved output ... need update 2012/12
vx_vols  = sqrt(sum(M0(1:3,1:3).^2));
vx_volt  = sqrt(sum(M1(1:3,1:3).^2));
interpol = any(vx_vols>vx_volt) + any(vx_vols/2>vx_volt);
if interpol
	yx = zeros(min([inf inf inf 3],size(yi)*2^interpol - 2^interpol + 1),'single');
	if interpol
		for i=1:3, yx(:,:,:,i) = single(interp3(yi(:,:,:,i),interpol,'cubic')); end % linear
	end
else
	yx=yi; % * regres/newres;
end

clear yi vx_vols vx_volt;

%% Saves the file
% Creates the output pathname
[Fpath,Fname] = spm_fileparts(nameOut);
nameOut  = fullfile(Fpath,'mri',['y_' Fname '_' num2str(numIteration) '.nii']);

% Prints a message about saving
fprintf('%s %d %s %s\n','Saving GS intermediate step',numIteration,'to a file:',nameOut);

% Creates and saves the NIFTI file
Yy2       = spm_diffeo('invdef',yx,odim,eye(4),M0);
N         = nifti;
N.dat     = file_array(nameOut,[odim(1:3),1,3],'float32',0,1,0);
N.mat     = M1;
N.mat0    = M1;
N.descrip = 'Deformation';
create(N);
N.dat(:,:,:,:,:) = reshape(Yy2,[odim,1,3]);
end
