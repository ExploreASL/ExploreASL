function [ EffectiveResolution ] = xASL_init_DefaultEffectiveResolution(ASLname, x)
%xASL_init_DefaultEffectiveResolution This function estimates what the ASL effective resolution
% should be based on simulations from Jan Petr

%% Admin

tIM         = xASL_io_ReadNifti(ASLname);
NativeRes   = tIM.hdr.pixdim(2:4);


% Determine estimated FWHM to convert native image resolution to
% effective resolution
if      strcmp(x.Sequence,'2D_EPI')
%         PAPER JAN PETR
%         EstimatedEffectiveResolution    = [3.2 3.2 7.4]; % Paper Jan Petr, 3.7 & 3.6 where X Y, which should be similar
%         BasedOnResolution               = [3 3 7];
%         Estimated_FWHM                  = EstimatedEffectiveResolution./BasedOnResolution;

        Estimated_FWHM                   = [1 1 1]; % KISS

elseif  strcmp(x.Sequence,'3D_spiral_WIP')
        % estimate
        Estimated_FWHM                   = [1.1158 1.1289 1.9525]; % average between 3D_spiral & 3D_GRASE


elseif  strcmp(x.Sequence,'3D_spiral')
%       PAPER JAN PETR
        EstimatedEffectiveResolution    = [4.3 4.4 10.1]; % Paper Jan Petr, 4.9 & 5.1 where X Y, which should be similar
        BasedOnResolution               = [3.8 3.8 4];
        Estimated_FWHM                  = EstimatedEffectiveResolution./BasedOnResolution;

        % The in-plane 2 translates the reconstruction resolution into
        % acquisition resolution

elseif  strcmp(x.Sequence,'3D_GRASE')
% %         Estimated_FWHM                  = ([4.5 4.5 6] + [2.6667 2.6667 2.3750]) ./2;
% %         % This is a rough educated guess, FWHM of 3D GRASE should be in between 3D spiral & 2D EPI.
%
        Estimated_FWHM                  = [1.1 1.1 1.38]; % KISS, Vidoretta, NeuroImage 2013
        % The in-plane 2 translates the reconstruction resolution into
        % acquisition resolution

% rough empirical guess!
%         EstimatedEffectiveResolution    = [4.5 4.5 6]; % Paper Jan Petr, 4.9 & 5.1 where X Y, which should be similar
%         BasedOnResolution               = [4   4   4];
%         Estimated_FWHM                  = EstimatedEffectiveResolution./BasedOnResolution;


end

if     ~isempty(regexp(x.Sequence,'3D_spiral')) && NativeRes(1)<4
        NativeRes(1:2)  = max(NativeRes(1:2),[3.75 3.75]);
        % to account for in-plane interpolation
elseif ~isempty(regexp(x.Sequence,'3D_GRASE')) && NativeRes(1)<3
        NativeRes(1:2)  = max(NativeRes(1:2),[3.8 3.8]);
end

EffectiveResolution     = NativeRes.*Estimated_FWHM;

fprintf('%s\n',[x.Sequence ' nifti has native resolution ' num2str(NativeRes(1)) ' ' num2str(NativeRes(2)) ' ' num2str(NativeRes(3)) ', assuming PSF ' num2str(Estimated_FWHM(1)) ' ' num2str(Estimated_FWHM(2)) ' ' num2str(Estimated_FWHM(3)) ' this gives estimated effective resolution ' num2str(EffectiveResolution(1)) ' ' num2str(EffectiveResolution(1)) ' ' num2str(EffectiveResolution(3))])



end
