function xASL_wrp_GSSaveIntermedTrans(y,idim,odim,rdim,M0,M1,R,M1t,M1r,fname,it)
% xASL_wrp_GSSaveIntermedTrans
%
% FORMAT: xASL_wrp_GSSaveIntermedTrans(y,idim,odim,rdim,M0,M1,R,M1t,M1r,fname,it)
%
% INPUT:
% y           - ...
% idim        - ...
% odim        - ...
% rdim        - ...
% M0          - ...
% M1          - ...
% R           - ...
% M1t         - ...
% M1r         - ...
% fname       - ...
% it          - ...
%
% OUTPUT: n/a
%
% DESCRIPTION: n/a
%
% EXAMPLE: n/a
% __________________________________
% Copyright 2015-2021 ExploreASL

% update for output resolution
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

% Create the filename
% VT0.fname == fname
[pth,nam] = spm_fileparts(fname);
fnameNew  = fullfile(pth,'mri',['y_' nam '_' num2str(it) '.nii']);
fprintf('%s %d %s %s\n','Saving intermediate step',it,'to file:',fnameNew);
Yy2       = spm_diffeo('invdef',yx,odim,eye(4),M0);
N         = nifti;
N.dat     = file_array(fnameNew,[odim(1:3),1,3],'float32',0,1,0);
N.mat     = M1;
N.mat0    = M1;
N.descrip = 'Deformation';
create(N);
N.dat(:,:,:,:,:) = reshape(Yy2,[odim,1,3]);
clear Yy2;

end
