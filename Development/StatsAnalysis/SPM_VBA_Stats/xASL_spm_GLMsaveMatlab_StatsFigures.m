function xASL_spm_GLMsaveMatlab_StatsFigures(diff_view_mean, H_ttestCONTRAST,x)
%xASL_spm_GLMsaveMatlab_StatsFigures Print Figures into files, using matlab figure stuff
    
    if  x.S.nSets>2
        iDiff   = 1; % print only the f-maps (subtractions don't make sense here, or could do the SD)
    else
        iDiff   = 2; % print only the significant subtraction maps
    end

%     for iDiff=1:length(diff_view_mean)

    H_ttest             = H_ttestCONTRAST;
    H_diff              = H_ttest .* diff_view_mean{iDiff};


    %% Determine window viewing/thresholding        
    hiThresh            = xASL_stat_MeanNan(nonzeros(H_diff)) + (4 * xASL_stat_StdNan(nonzeros(H_diff)));
    hiThresh            = ceil(hiThresh/2)*2;
    loThresh            = xASL_stat_MeanNan(nonzeros(H_diff)) - (4 * xASL_stat_StdNan(nonzeros(H_diff)));
    loThresh            = floor(loThresh/2)*2;        

    if      x.S.nSets==1 && ~S.RegressionCOVAR % 1-sample t-test
            ThreshPlus              = hiThresh;
            ThreshMin               = loThresh-(hiThresh-loThresh);
    elseif  x.S.RegressionCOVAR % regression
            ThreshPlus              = hiThresh;
            ThreshMin               = loThresh;


    else % Obtain reasonable clipping values for viewing (window level)
            % Window level should be same for increase & decrease
            ThreshStat              = max(abs([hiThresh loThresh]));
            ThreshPlus              = ThreshStat;
            ThreshMin               = -ThreshStat;
    end

    if      isfield(x,'StatWindowLevel')
            ThreshStat              = x.StatWindowLevel;
            ThreshPlus              = ThreshStat;
            ThreshMin               = -ThreshStat;
    end


    
    %% Clip data
    H_diff(H_diff<ThreshMin)        = ThreshMin;
    H_diff(H_diff>ThreshPlus)       = ThreshPlus;
    clear diff_view hiThresh loThresh ThreshStat

    % Clip at window level
    H_diff(isnan(H_diff))       = 0;

    close all

    %% Create colormap
    fig = figure('Visible','off');
    subplot(1,1,1);

    joined_colormap     = xASL_vis_joinColormap(1,x.S.cool,x.S.hot);
    joinedTickLabels    = [ThreshMin ThreshPlus]; % = get(handleBar,'TickLabels'); % this doesn't work
    
    if  xASL_stat_SumNan(H_ttest(:))>0 % if there are significant differences, otherwise skip to below, to show background only
        imshow(  H_diff,[ThreshMin ThreshPlus],'border','tight');
		
		% Get the handle of the axes created by imshow
		handleAxes          = get(get(0,'CurrentFigure'),'CurrentAxes');
		
		% change the colormap of the axes
		colormap(handleAxes,joined_colormap);
        handleBar           = colorbar();
		
        % Convert significance map to joined_colormap
        % First split H_diff into negative & positive difference map
        % Colormap runs from 0:255 instead of 1:256,
        % So 127 instead of 128 represents 0.

        NewMatrix         = xASL_im_RescaleColorbarStats( H_diff);
        NewMatrixClr      = ind2rgb(double(NewMatrix),joined_colormap);
    else
        NewMatrixClr      = repmat(zeros(size(H_diff)),[1 1 3]);
    end

    
    
    %% Recreate background image
    [ x ]                                  = xASL_init_PopulationSettings( x );
    
    MaskBackground                         = repmat(H_diff==0,[1 1 3]); % without any significant cluster, this shows background only
    MaskBackground(isnan(MaskBackground))  = 0;
    MaskBackground                         = logical(MaskBackground);
    NewMatrixClr(MaskBackground)           = x.background_view_clr(MaskBackground);

    warning('off','images:initSize:adjustingMag'); % warning about scaling if an image doesnt fit screen, disable

    % Show new composed map
    imshow(NewMatrixClr,'border','tight');
    axesBar = xASL_vis_Colorbar(joined_colormap,joinedTickLabels);
    
    %% Text stuff
    printTitle  = x.S.printTitleORI;

    if          iDiff==1 % stat-map
                printTitle          = [printTitle ' t_f-map'];
    end

    if  isfield(x.S,'PerSetName') % Add previous set-selection if any ->
%            kept out because already in pathname, and would make filename too
%            long
            printTitle      = [printTitle ' ' x.S.PerSetName];
    end

    if  exist('CoVariate','var')
        if      x.S.RegressionCOVAR
                printTitle      = [printTitle ' predictor '  CoVariate.Name];
        else
                printTitle      = [printTitle ' co-variate ' CoVariate.Name];
        end
    end

    title( printTitle,'interpreter','none','HorizontalAlignment','right');
    xlabel( x.S.PrintTitleStats,'interpreter','none','HorizontalAlignment','right'  );

    if  x.S.GlobalNormalization
        x.S.unit  = ['%' x.S.unit];
    end

	if      x.S.nSets>2
		axes(axesBar);ylabel('Coefficient of Variation (100x(SD/mean))');
	elseif  x.S.RegressionCOVAR
		axes(axesBar);ylabel(['Beta '  x.S.OutputID ' (' x.S.unit ')']);
	elseif  x.S.nSets==1
		axes(axesBar);ylabel(['Mean '  x.S.OutputID ' (' x.S.unit ')']);
	else
		axes(axesBar);ylabel(['delta ' x.S.OutputID ' (' x.S.unit ')']);
	end

    %% Print results
    % Tried several options here.
    % Creating a pdf will still not be optimal for illustrator, and
    % takes about 6 seconds. Creating an EPS takes little time,
    % creating a TIFF takes little time as well.
    % The EPS contains very high resolution images, the TIFF resolution
    % can be adjusted to taste, 300 DPI is a good tradeoff between
    % quality & speed/datasize for larger figures.

    SaveFile           = fullfile(x.S.StatsDir,[printTitle '.tiff']);
    print('-dtiff','-r300',SaveFile); % saveas( fig , SaveFile,'jpg');
    SaveFile           = fullfile(x.S.StatsDir,[printTitle '.eps']);
    saveas( fig ,SaveFile,'epsc');
    close all

   
end




%% ===========================================================================================
%% ===========================================================================================

function [ NewMatrix ] = xASL_im_RescaleColorbarStats( OriMatrix)
%xASL_im_RescaleColorbarStats Linearly rescales input matrix to colorbar output matrix,
% where [  0:1:127] represents negative changes (<0)
% and   [128:1:255] represents positive changes (>0)
% 128 is assumed to be 0 on the colorbar

% Split map into negative & positive maps
NegativeMatrix    = OriMatrix.*(OriMatrix<0);
PositiveMatrix    = OriMatrix.*(OriMatrix>0);

% Linearly rescale both maps
negResc     = round(xASL_vis_RescaleLinear(NegativeMatrix,  1,128,0));
posResc     = round(xASL_vis_RescaleLinear(PositiveMatrix,129,256,0));

% Combine maps
NewMatrix                   = zeros(size(OriMatrix,1),size(OriMatrix,2));
NewMatrix(:)                = 128;
NewMatrix(negResc<128)      = negResc(negResc<128);
NewMatrix(posResc>129)      = posResc(posResc>129);

end


function [ NewMatrix ] = xASL_vis_RescaleLinear(OriMatrix,NewMin,NexMax,NonZerosOption)
%rescale Linearly rescales input matrix to output matrix,
% applying a new minimum and new maximum.

% % Test values
% % old_matrix      =round(rand(5,5).*10);
% % new_minimum     =2;
% % new_maximum     =6;
% % 

OriMatrix  = double( OriMatrix );

% If matrix contains zeros only, reinforce nonzerosoption==0
if  max(OriMatrix(:))==0
    NonZerosOption=0;
end

% define minimum old matrix
if      NonZerosOption==1
        OriMin         = min(nonzeros(OriMatrix(:)));
else
        OriMin         = min(OriMatrix(:));
end


OriMax                 = max(OriMatrix(:));
OriRange               = OriMax-OriMin;
NewRange               = NexMax-NewMin;

NewMatrix              = NewMin+NewRange.* (OriMatrix-OriMin)./OriRange;


end
