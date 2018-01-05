function visualizeSignalCoefficients2Dmap(ax, recCoef, WFScoef, NScoef)

powRec = (abs(recCoef).^2)/2;
powWFS = (abs(WFScoef).^2)/2;
powNS = (abs(NScoef).^2)/2;

powRecDB = 10*log10(powRec);
powWFSDB = 10*log10(powWFS);
powNSDB = 10*log10(powNS);

scatRec = findobj(ax, 'Tag', 'receiver');
scatWFS = findobj(ax, 'Tag', 'loudspeakers');
scatNS = findobj(ax, 'Tag', 'source');

% % Set RGB color
% recColorMap = colormap('gray');
% indices = scaled2indexedColors(size(recColorMap, 1), [-100 20], powRecDB);
% scatRec.CData = recColorMap(indices, :);
% 
% WFSColorMap = colormap('gray');
% indices = scaled2indexedColors(size(WFSColorMap, 1), [-100 20], powWFSDB);
% scatWFS.CData = WFSColorMap(indices, :);

% Set color and 
ax.CLim = [-100, 20];
scatRec.CData = powRecDB;
scatWFS.CData = powWFSDB;
scatNS.CData = powNSDB;

end