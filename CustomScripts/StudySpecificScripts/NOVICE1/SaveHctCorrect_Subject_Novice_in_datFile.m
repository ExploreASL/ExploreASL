opefunction SaveHctCorrect_SubjectSession_SLEEP( x)
%SaveHctCorrect Corrects CBF-values for individual hematocrit or group hematocrit
%
% By HJMM Mutsaerts, ExploreASL 2016
% NB: this script assumes a single hct value for each subject (for all sessions)
%
%

    load( fullfile( x.D.ROOT, 'Mean_Hct_Session_Cohort' ) );
    HCTALL  = HCT;
    clear HCT




    %% 2    Create CBF maps corrected for individual Hct values

    clear ASLnew

    for iSubject=1:x.nSubjects
        iSession    = 1;
            WasFound        = 0;
            SubjectSessN    = iSubject;
            HCTvalue        = HCTALL(SubjectSessN,1);


                        %% convert to numeric if not numeric
                        if      HCTvalue>1 && HCTvalue<100 % get percentage
                                HCTvalue  = HCTvalue / 100;
                        end


                        %% Create corrected ASL maps

                        FNsliceGradient     = fullfile( x.D.PopDir, ['DARTEL_slice_gradient_' x.SUBJECTS{iSubject} '_' x.SESSIONS{iSession} '.nii']);

                        SliceGradient       = xASL_io_ReadNifti( FNsliceGradient );
                        SliceGradient       = SliceGradient.dat(:,:,:);

                        SliceGradient       = xASL_im_ndnanfilter(SliceGradient,'gauss',[1.885 1.885 1.885],1);
                        tempIm              = squeeze(ASL.Data.data(SubjectSessN,:,:,:) );

%                         % Overlap check, since these
%                         slice_maps were from previous run.
%                         But they're still good.
%                         tempIm(isnan(tempIm))     = 0;
%                         dip_image(tempIm+(SliceGradient.*5))
%                         tempImOLD   = tempIm;

                        %% Create PLD gradient
                        if      isnumeric( x.Q.SliceReadoutTime )
                                PLDslicereadout     = x.Q.SliceReadoutTime;
                        end
                        SliceGradient  = x.Q.Initial_PLD + ((SliceGradient-1) .* PLDslicereadout); % effective PLD

                        %% Calculate blood T1 with hematocrit (Hales et al., 2016)
                        T1aNew      = calc_blood_t1(HCTvalue, 0.97, 3) .* 1000; % with Y=0.97 and for B0=3T

                        lab_eff     = 0.83*0.85;

                        QntFactorOld    = exp(SliceGradient./x.Q.BloodT1) / (2.*lab_eff.*x.Q.BloodT1 .* (1- exp(-x.Q.LabelingDuration./x.Q.BloodT1)) );
                        QntFactorNew    = exp(SliceGradient./T1aNew ) / (2.*lab_eff.*T1aNew  .* (1- exp(-x.Q.LabelingDuration./T1aNew )) );

                        tempIm          = (tempIm .* QntFactorNew) ./ QntFactorOld;
                        ASLnew(SubjectSessN,:,:,:)  = tempIm;
%                         xASL_io_SaveNifti( FileNameCBF, FileNameCBFCorrInd, tempIm );

                        clear T1aNew QntFactorOld QntFactorNew tempIm SliceGradient FileNameCBFCorrInd FNsliceGradient FileNameCBF
                        clear HCTvalue HCT HB B0 FE Y R1BL deltaT1BLOOD lab_eff SubjectSessN PLDslicereadout

                        WasFound = WasFound+1;
    end



    %

        % Notify missing hematocrit values
        if      WasFound==0
                fprintf('%s\n',[x.SUBJECTS{iSubject} ' was not found!!!']);
        elseif  WasFound>1
                fprintf('%s\n',['Duplicate hct values found for ' x.SUBJECTS{iSubject} '!!!']);
        end
        clear WasFound

        %% Save them


        fprintf('\n');
        LoadFile  = fullfile( x.D.PopDir, 'ASL_HctCohort.dat');
        memcheck_ASL;


        if ~exist( LoadFile ,'file' )
            % Temporary load all ASL files, to create single 5D matrix
            % CAVE: should be explicit file definition!
            % CAVE: memory mapping fails if exclusions are not explicitly defined! All inclusions should be available!

            % Pre-allocate memory
            ASL_HctCohort   = single(ASLnew);

            % clip for viewing
            ASL_HctCohort(ASL_HctCohort<0)      = 0;

            fileID = fopen(LoadFile,'w');
            fwrite( fileID,ASL_HctCohort,'single');
            fclose(fileID);
        end


end
