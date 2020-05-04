function x    = xASL_wrp_VisualCheckCollective_ASL(x)
%xASL_wrp_VisualCheckCollective_ASL Runs a collection of visual QC
%functions




%% Get visualization settings
% Parameters for creating visual QC Figures:
% CBF, CBF with overlay c2T1, CBF with overlay c2T1
% MeanControl SD
% M0 NoSmoothM0 NoSmoothM0 with overlay c1T1
% TT TT with overlay c2T1

x               = xASL_adm_ResetVisualizationSlices(x);

T.ImIn          = {x.P.Pop_Path_qCBF  x.P.Pop_Path_SD {x.P.Pop_Path_qCBF x.P.Pop_Path_PV_pWM} x.P.Pop_Path_SNR};
T.ImIn(5:8)     = {x.P.Pop_Path_mean_control x.P.Pop_Path_noSmooth_M0 {x.P.Pop_Path_noSmooth_M0 x.P.Pop_Path_PV_pGM} x.P.Pop_Path_M0};
T.ImIn(9:10)    = {x.P.Pop_Path_TT  {x.P.Pop_Path_TT x.P.Pop_Path_PV_pWM}};

T.DirOut        = {x.D.ASLCheckDir x.D.SNRdir      x.D.ASLCheckDir       x.D.SNRdir};
T.DirOut(5:8)   = {x.D.RawDir      x.D.M0CheckDir  x.D.M0regASLdir       x.D.M0CheckDir};
T.DirOut(9:11)  = {x.D.TTCheckDir  x.D.TTCheckDir  x.D.ASLCheckDir};

T.IntScale(2)   = {[1 1]};
T.IntScale{8}   = [0.75 0.65];

T.ColorMapIs{10}= x.S.jet256;
T.ColorMapIs{11}= {x.S.jet256};

T.NameExt       = {[] [] 'Reg_pWM' []};
T.NameExt( 5:8) = {[] [] 'Reg_pGM' []};
T.NameExt(9:11) = {[] 'Reg_pWM' 'AnalysisMask'};

% Fill missing cells
Pars = {'ImIn' 'DirOut' 'ClipZero' 'IntScale' 'NameExt' 'ColorMapIs'}; % default pars
for iM=1:length(T.ImIn)
    for iP=1:length(Pars)
        if ~isfield(T,Pars{iP})
            T.(Pars{iP}) = [];
        elseif length(T.(Pars{iP}))<iM
            T.(Pars{iP}){iM} = [];
        end
    end
end

%% Take only NIfTIs that exist
for iL=1:length(T.ImIn)
    if  ischar(T.ImIn{iL})
        ExistInd(iL) = logical(xASL_exist(T.ImIn{iL},'file'));
    else
        ExistInd(iL) = min(cellfun(@(x) logical(xASL_exist(char(x),'file')),T.ImIn{iL}));
    end
end
        
for iP=1:length(Pars)
    T.(Pars{iP}) = T.(Pars{iP})(ExistInd);
end

%% Clone each row into a transversal & coronal row
nIms = length(T.(Pars{1}));
nRows = ceil( nIms/4);

for iN=1:nRows
    clear T2
    ImsI                    = (iN-1)*4+1:min(nIms,iN*4);
    nImsRow                 = length(ImsI);
    nRow1                   = 1:nImsRow;
    nRow2                   = nImsRow+1:2*nImsRow;
    T2.TraSlices(nRow1)     = {[]};
    T2.TraSlices(nRow2)     = {'n/a'};
    T2.CorSlices(nRow1)     = {'n/a'};
    T2.CorSlices(nRow2)     = {[]};

    T2.NameExt(nRow1)       = cellfun(@(x) ['Tra_' x], T.NameExt(ImsI), 'UniformOutput',false);
    T2.NameExt(nRow2)       = cellfun(@(x) ['Cor_' x], T.NameExt(ImsI), 'UniformOutput',false);
    T2.ImIn                 = [T.ImIn(ImsI) T.ImIn(ImsI)];
    T2.DirOut               = [T.DirOut(ImsI) T.DirOut(ImsI)];
    T2.IntScale             = [T.IntScale(ImsI) T.IntScale(ImsI)];
    T2.ColorMapIs           = [T.ColorMapIs(ImsI) T.ColorMapIs(ImsI)];
    T2.ModuleName           = 'ASL';

%%  Perform the visualization
% Perhaps at the end of the row we need to generate empty images, as transversal & coronal have different sizes
% they don't concatenate well horizontally, need to be concatenated vertically

    fprintf('%s','Printing images...  ');
    for iM=1:length(T2.ImIn)
        xASL_TrackProgress(iM,length(T2.ImIn)*nRows);

        % Manage slices to show
        % Sagittal
        x.S.SagSlices   = []; % show no sagittal slices
        % Transversal
        if      isempty(T2.TraSlices{iM})
                x.S.TraSlices   = x.S.slicesLarge;
        elseif  strcmp(T2.TraSlices{iM},'n/a')
                x.S.TraSlices   = [];
        else
                warning('Wrong slice choice');
        end
        % Coronal
        if      isempty(T2.CorSlices{iM})
                x.S.CorSlices   = x.S.slicesLarge+7;
        elseif  strcmp(T2.CorSlices{iM},'n/a')
                x.S.CorSlices   = [];
        else
                warning('Wrong slice choice');
        end    

        % Create the image
        T2.IM = xASL_im_CreateVisualFig( x, T2.ImIn{iM}, T2.DirOut{iM}, T2.IntScale{iM}, T2.NameExt{iM}, T2.ColorMapIs{iM});
        % add single slice to QC collection
        if sum(~isnan(T2.IM(:)))>0 % if image is not empty
            x = xASL_im_AddIM2QC(x,T2);
        end
    end
end

fprintf('\n');
   
end
