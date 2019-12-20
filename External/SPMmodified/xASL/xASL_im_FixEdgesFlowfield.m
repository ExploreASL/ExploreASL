function xASL_im_FixEdgesFlowfield(InputPath)
%xASL_im_FixEdgesFlowfield Iteratively smooths the flowfields to its edges,
% filling any artifacts created from interpolating NaNs, at combining flowifleds

    fprintf('%s', 'Cleaning up flowfield edges:  0%');

    InputIM = xASL_io_Nifti2Im(InputPath);

    for iM=1:3
        xASL_TrackProgress(iM-1, 3);
        IM = InputIM(:,:,:,:,iM);
        while sum(isnan(IM(:)))>0
               IM = xASL_im_ndnanfilter(IM, 'gauss', [8 8 8], 2);
        end
        InputIM(:,:,:,:,iM) = IM;
        xASL_TrackProgress(iM, 3);
    end
    fprintf('\n');

    xASL_io_SaveNifti(InputPath, InputPath, InputIM, [], 0);

end
