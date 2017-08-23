% Import image
imWFSArray = imread('WFSarrayScheme.png', 'png', 'BackgroundColor', [1 1 1]);

xLim = [viewBox(1), viewBox(1) + viewBox(3)];
yLim = [viewBox(2), viewBox(2) + viewBox(4)];

ax = axes(figure);
image(ax, 'XData', xLim, 'YData', yLim, 'CData', imWFSArray);
ax.DataAspectRatio = [1, 1, 1];

