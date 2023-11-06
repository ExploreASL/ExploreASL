function xASL_io_PairwiseSubtraction(InputFile,outputPath,do_mask,switch_sign)
% xASL_io_PairwiseSubtraction Subtracts controls from labels and takes mean
% Creates new perfusion-weighted delta_M file, prefaced with 's'
% Converts into single precision floating point values (32 bit), removes scale slope.
% Only runs if ['s' input_file_ASL] doesn't exist
% Remember to consider motion correction/ SPM realign (isotropically)
% Alternative to this function is robust fit (Camille Maumet)
%
% FORMAT:       xASL_io_PairwiseSubtraction(InputFile,outputPath,do_mask,switch_sign)
% 
% INPUT:        do_mask       = masks perfusion-weighted image based on voxels present in all temporal volumes (default=0)
%               switch_sign   = script checks whether odd or even are largest (control),
%                               in order to calculate (control-label).
%                               if this nevertheless goes wrong, switch_sign=1 will invert this. (default=0)
%
%               Script dependencies: SPM12
%               Default input values.
%               If dim(4)==1, it will not process anything but rather
%               create a copy of the original nifti.
%
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Subtracts controls from labels and takes mean.
%               Creates new perfusion-weighted delta_M file, prefaced with 's'.
%               Converts into single precision floating point values (32 bit), removes scale slope.
%               Only runs if ['s' input_file_ASL] doesn't exist.
%               Remember to consider motion correction/ SPM realign (isotropically).
%               Alternative to this function is robust fit (Camille Maumet).
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL


if ~exist('do_mask','var'); do_mask             = 0; end
if ~exist('switch_sign','var'); switch_sign     = 0; end

% Redefine output_file, to make sure that it doesn't contain path
[path file ext]     = xASL_fileparts(outputPath);
OutputFile          = [file ext]; % If there is no path, output_file will remain same
ControlFile         = 'mean_control.nii';
OutputPath2         = fullfile(path,ControlFile);

if ~xASL_exist(InputFile,'file')
        warning([InputFile ' did not exist']);
else

    xASL_delete( outputPath);
    xASL_delete(OutputPath2);

    tIM             = xASL_io_Nifti2Im(InputFile);

    if  length(size(tIM))<4 % if it concerns a single frame, simply copy it
        %xASL_Copy(InputFile, outputPath,1);
		tIM(isnan(tIM)) = 0;
		xASL_io_SaveNifti(InputFile,outputPath,tIM,16,0);
        fprintf('%s\n','mean_control.nii not created because of single frame, no timeseries');
    else
        maskvalue       = 0.2.*max(tIM(:));


        % Check order control-label frames
        [control_im, label_im, OrderContLabl] = xASL_quant_GetControlLabelOrder(tIM);

        % Check equality of n frames control & label
        if  size(control_im,4)~=size(label_im,4)
            error('Control and label image have unequal n frames!');
        end

        % switch_sign & OrderContLabl should be -1 or 1
        if ~switch_sign
            switch_sign     = -1; % switch_sign==1 is change order
        end
        if ~OrderContLabl
            OrderContLabl   = -1; % OrderContLabl==1 is change order
        end

        % if both 0 or 1, then leave order ((-1 * -1) > 0)
        % if one 0 other 1, then change order ((-1 * 1) < 0)

        % default is control label control label

        if  (OrderContLabl.*switch_sign)<0 % change order if necessary
             largest 		= 2;
        else largest 		= 1;
        end

        PWI_expr            = '0';
        raw_expr            = '0';
        mask                = '0';
        Incl                = 0;
        excluded            = 0;
        excluded_pairs      = '';


       % Create equation, with all volumes

        for iI=1:size(tIM,4)/2
           if  (OrderContLabl.*switch_sign)<0 % change order
              PWI_expr = [PWI_expr '-i' num2str(iI*2-1) '+i' num2str(iI*2)];
              raw_expr = [raw_expr '+i' num2str(iI*2)];
           else
              PWI_expr = [PWI_expr '+i' num2str(iI*2-1) '-i' num2str(iI*2)];
              raw_expr = [raw_expr '+i' num2str(iI*2-1)];
           end
           mask         = [mask '+(i' num2str(iI*2-1) '>' num2str(maskvalue) ')+(i' num2str(iI*2) '>' num2str(maskvalue) ')'];
           Incl         = Incl+1;
        end

        PWI_expr    = ['((' PWI_expr ')/' num2str(Incl) ')'];
        raw_expr    = ['((' raw_expr ')/' num2str(Incl) ')'];

        if  do_mask
            PWI_expr    = [PWI_expr ' .* ((' mask ')>=' num2str(Incl) ')'];
            raw_expr    = [raw_expr ' .* ((' mask ')>=' num2str(Incl) ')'];
        end


        % Other imcalc parameters
        matlabbatch{1}.spm.util.imcalc.outdir           = {path};
        matlabbatch{1}.spm.util.imcalc.options.dmtx     = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask     = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp   = 4;
        matlabbatch{1}.spm.util.imcalc.options.dtype    = 16;
        for ii=1:size(tIM,4); matlabbatch{1}.spm.util.imcalc.input{ii,1} = [InputFile ',' num2str(ii)]; end  % Input all volumes


        %% Run job
        matlabbatch{1}.spm.util.imcalc.output           = OutputFile;
        matlabbatch{1}.spm.util.imcalc.expression       = PWI_expr;
        fprintf('%s\n',['Computing ' OutputFile]);
        spm_jobman('run', matlabbatch);

        matlabbatch{1}.spm.util.imcalc.output           = ControlFile;
        matlabbatch{1}.spm.util.imcalc.expression       = raw_expr;
        fprintf('%s\n',['Computing ' ControlFile]);
        spm_jobman('run', matlabbatch);
    end

end
