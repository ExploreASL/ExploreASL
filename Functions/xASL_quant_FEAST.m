function xASL_quant_FEAST(x)
%xASL_quant_FEAST Computation FEAST-based transit times (uses images that were not vasculary treated)
% if there is a crushed & non-crushed scan, then transit times will be
% computed by division of these scans, provided sessions are exactly named
% as defined below

if x.nSessions>1 && isfield(x,'session') && isfield(x.session,'options') && strcmp(x.session.options{1},'crushed') && strcmp(x.session.options{2},'non-crushed')
    if strcmp(x.SESSIONDIR(length(x.SUBJECTDIR)+2:end),'ASL_2') % Computation is performed if CurrentSession=session 2
        for iSession=1:2
            PathCBF{iSession} = fullfile(x.D.PopDir, ['qCBF_' x.P.SubjectID '_' x.SESSIONS{iSession} '.nii']);
        end
        if ~xASL_exist(PathCBF{1},'file') || ~xASL_exist(PathCBF{2},'file')
            return; % skip if files dont exist (yet)
        end

        fprintf('%s\n','Saving TT nifti');

        xASL_adm_CreateDir(x.D.TTCheckDir);

        %% load data

        for iSession=1:2
            % Load data
            CBF{iSession} = xASL_io_Nifti2Im(PathCBF{iSession});
            Gradient{iSession} = xASL_io_Nifti2Im(fullfile(x.D.PopDir, ['SliceGradient_extrapolated_' x.P.SubjectID '_' x.SESSIONS{iSession} '.nii']));
            % correct different PLD scales
            PLD{iSession} = x.Q.Initial_PLD + ((Gradient{iSession}-1) .* x.Q.SliceReadoutTime);
            Gradient{iSession} = exp(PLD{iSession}./x.Q.BloodT1) / (2.*x.Q.LabelingEfficiency.*x.Q.BloodT1 .* (1- exp(-x.Q.LabelingDuration./x.Q.BloodT1)) );
            CBF{iSession} = CBF{iSession}./Gradient{iSession};
        end
        % Average different PLD scales
        PLD_combined= (PLD{1}+PLD{2})./2;

        for iSession=1:2
            % Masking
            CBF{iSession}(~x.skull) = NaN; % this is not masked for WM or GM, we could do that
            % smooth maps, ignoring NaNs
            % CAVE: NaNs are interpolated with data, hence the maps should be masked later!
            
            CBF{iSession} = xASL_im_ndnanfilter(CBF{iSession},'gauss',[8 8 8],0);
        end
        FEAST_ratio = CBF{1}./CBF{2}; % crushed/non-crushed
        FEAST_ratio(FEAST_ratio>1) = 1; % clip @ 1 (== ATT = PLD)
        FEAST_ratio(FEAST_ratio<0) = 0; % clip @ 0 (== infinite ATT)

        % Calculate 
        qnt_PLDdecay = exp(-PLD_combined/x.Q.BloodT1);
        qnt_combidecay = exp( (-x.Q.LabelingDuration - PLD_combined) / x.Q.BloodT1);
        TT = -x.Q.BloodT1 .* reallog( FEAST_ratio .* (qnt_PLDdecay  - qnt_combidecay ) + qnt_combidecay );

        xASL_io_SaveNifti(PathCBF{1}, fullfile(x.D.PopDir, ['TT_'  x.P.SubjectID '.nii']), TT, 32);
    end
end


end

