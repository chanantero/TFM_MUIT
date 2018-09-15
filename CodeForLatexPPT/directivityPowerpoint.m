%% Monopole and Dipole for Powerpoint presentation

%% Simple image

backgroundFileName = [imagesPath, 'KirchhoffTheoScheme_ppt_SurfaceDif.png'];
[A,map,transparency] = imread(backgroundFileName);
[imHeight, imWidth, numColors] = size(A);

% axBack = axes(figure);
% image(axBack, A, 'AlphaData', transparency)

refImage = [417, 302]; % [x, y] Píxel coordinates where the origin of coordinates will be
imScale = 200; % Number of pixels that represent a unit of length

XLim = [-1, 1];
YLim = [-1, 1];
imXLim = [-refImage(1), imWidth - 1 - refImage(1)]/imScale;
imYLim = [-refImage(2), imHeight - 1 - refImage(2)]/imScale;
axBack = axes(figure, 'NextPlot', 'Add');
axBack.DataAspectRatio = [1 1 1];
image(axBack, imXLim, imYLim, A, 'AlphaData', transparency)

% Dipole
alpha = linspace(0, 2*pi, 100);
dipoleMag = abs(cos(alpha));
dipoleColor = [237, 125, 49]/255;
[xDip, yDip] = pol2cart(alpha, dipoleMag);

axDip = copyobj(axBack, figure);
plot(axDip, xDip, yDip, 'Color', dipoleColor, 'LineWidth', 3)
axDip.Visible = 'off';
scatter(axDip, 0, 0, 50, dipoleColor, 'filled')

saveas(axDip.Parent, [imagesPath, 'KirchhoffTheoScheme_ppt_Dipole.svg'], 'svg')

% Monopole
monopoleMag = ones(size(alpha));
monopoleColor = [68, 114, 196]/255;
[xMon, yMon] = pol2cart(alpha, monopoleMag);

axMon = copyobj(axBack, figure);
plot(axMon, xMon, yMon, 'Color', monopoleColor, 'LineWidth', 3)
axMon.Visible = 'off';
scatter(axMon, 0, 0, 50, monopoleColor, 'filled')

saveas(axMon.Parent, [imagesPath, 'KirchhoffTheoScheme_ppt_Monopole.svg'], 'svg')


% Both
axBoth = copyobj(axBack, figure);
plot(axBoth, xMon, yMon, 'Color', monopoleColor, 'LineWidth', 3)
plot(axBoth, xDip, yDip, 'Color', dipoleColor, 'LineWidth', 3)
axBoth.Visible = 'off';

saveas(axBoth.Parent, [imagesPath, 'KirchhoffTheoScheme_ppt_Both.svg'], 'svg')


