function axesBar = xASL_io_colorbar(ColorMap,TickLabels)

% Display a color bar, but replace it with a colorbar with given 'ticks' and colormap

% Display colorbar and get its handle
handleBar = colorbar();
propsBar = get(handleBar);

% Create new axes
handleAx = axes();

% Display the colormap in the new axes
hold('all')	
imagesc(ind2rgb((1:256)',ColorMap),'Parent',handleAx);
set(handleAx,'YLim',[1 256]);

% Adjust to the size of the original colorbar, with correct ticks and location.
set(handleAx,'Position',propsBar.Position);

YTicks          = round(1:25:256);
TickSpread      = max(TickLabels)-min(TickLabels);
TickSpacing     = TickSpread/(length(YTicks)-1);
NewLabels       = [min(TickLabels):TickSpacing:max(TickLabels)];

set(handleAx,'YTickLabel',NewLabels);
set(handleAx,'YTick',YTicks);
set(handleAx,'XTick',[]);
set(handleAx,'YAxisLocation','right')

% Remove the original colorbar.
set(handleBar,'Visible','off');

hold('off')

axesBar = handleAx;
return