%% Draw WFS theory schemes
% pathSetUp;
% imagesPath = 'C:/Users/Rubén/Google Drive/Telecomunicación/Máster 2º Curso 2015-2016/TFM MUIT/Documentos/TFM/Img/';

%% Kirchhoff dimensionality reduction

ax = axes(figure, 'NextPlot', 'Add');
xlabel('x')

% Create WFS plane
x = [0 1 1 0 0];
y = [0 0 1 1 0];
z = zeros(size(x));
plot3(ax, x, y, z, 'k')

% Create section of the surface S as an ellipse
x1 = 0.3; x2 = 0.9; y1 = 0.9; y2 = 0.3; % Coordinates of the two extremes of the major axis of the ellipse
e = 0.8; % Excentricity
a = 1/2*sqrt((x2-x1)^2+(y2-y1)^2);
b = a*sqrt(1-e^2);
t = linspace(0,2*pi, 100);
X = a*cos(t);
Y = b*sin(t);
w = atan2(y2-y1,x2-x1);
xSect = (x1+x2)/2 + X*cos(w) - Y*sin(w);
ySect = (y1+y2)/2 + X*sin(w) + Y*cos(w);
zSect = zeros(size(xSect));
plot3(ax, xSect, ySect, zSect, 'k')

% Extend the tube
h = 10;
s = surf(ax, [xSect; xSect], [ySect; ySect], [zSect + h; zSect - h]);
s.LineStyle = 'none';

% Lighting properties
s.FaceLighting = 'gouraud';
s.BackFaceLighting = 'lit';
s.FaceNormals = -s.FaceNormals;
l = light(ax, 'Position', [-1 0 0], 'Style', 'local');

% Create primary source
PSpos = [0.1 0.1 0];
ps = scatter3(ax, PSpos(1), PSpos(2), PSpos(3), 50, [0 0 0], 'filled');

% Create receiving point
RecPos = [0.7 0.5 0];
recP = scatter3(ax, RecPos(1), RecPos(2), RecPos(3), 50, [0 0 0], 'filled');

% Create normal arrow to S in section at a random point
ind = 60;
SSpos = [xSect(ind), ySect(ind), zSect(ind)]; % Secondary source position
tangentDir = atan2d(ySect(ind+1) - ySect(ind), xSect(ind+1) - xSect(ind));
normalDir = tangentDir + 90;
longArrow = 0.15;
arrow([SSpos(1) SSpos(2) 0], [SSpos(1) + longArrow*cosd(normalDir), SSpos(2) + longArrow*sind(normalDir), 0],...
    'CrossDir', [1 0 0], 'NormalDir', [0 0 1])

% Create r0 vector
plot3(ax, [PSpos(1), SSpos(1)], [PSpos(2), SSpos(2)], [PSpos(3), SSpos(3)], 'k');
drawbrace([SSpos(1) SSpos(2)], [PSpos(1) PSpos(2)], 10); % Add curly brace
dir = [SSpos(1) - PSpos(1), SSpos(2) - PSpos(2), SSpos(3) - PSpos(3)]; dir = dir/norm(dir);
plot3(ax, [SSpos(1), SSpos(1) + dir(1)*longArrow], [SSpos(2), SSpos(2) + dir(2)*longArrow], [SSpos(3), SSpos(3) + dir(3)*longArrow], 'k--');

% Create Deltar0 vector
plot3(ax, [RecPos(1), SSpos(1)], [RecPos(2), SSpos(2)], [RecPos(3), SSpos(3)], 'k');
drawbrace([RecPos(1) RecPos(2)], [SSpos(1) SSpos(2)], 10); % Add curly brace

% Draw arc for the incidence angle
center = SSpos;
radius = 0.1;
alpha1 = atan2d(SSpos(2) - PSpos(2), SSpos(1) - PSpos(1));
alpha2 = normalDir;
alpha = linspace(alpha1, alpha2, 10);
xArc1 = center(1) + radius*cosd(alpha);
yArc1 = center(2) + radius*sind(alpha);
zArc1 = zeros(size(xArc1));
posMidArc1 = [center(1) + radius*cosd((alpha1 + alpha2)/2), center(2) + radius*sind((alpha1 + alpha2)/2)];
plot3(ax, xArc1, yArc1, zArc1, 'k:')

% Draw arc for the propagation angle from the secondary source to the receiving point
center = SSpos;
radius = 0.08;
alpha1 = atan2d(RecPos(2) - SSpos(2), RecPos(1) - SSpos(1));
alpha2 = normalDir;
alpha = linspace(alpha1, alpha2, 10);
xArc2 = center(1) + radius*cosd(alpha);
yArc2 = center(2) + radius*sind(alpha);
zArc2 = zeros(size(xArc1));
posMidArc2 = [center(1) + radius*cosd((alpha1 + alpha2)/2), center(2) + radius*sind((alpha1 + alpha2)/2)];
plot3(ax, xArc2, yArc2, zArc2, 'k:')

ax.Visible = 'off';

% Add latex labels
t = gobjects(1);
offsetX = 0.02; offsetY = -0.02;
t(1) = text(ax, PSpos(1) + offsetX, PSpos(2) + offsetY, '$\PosTheo[primarySource]$');
t(2) = text(ax, SSpos(1) + offsetX, SSpos(2) + offsetY, '$\PosTheo[section]$');
t(3) = text(ax, RecPos(1) + offsetX, RecPos(2) + offsetY, '$\PosTheo$');
t(4) = text(ax, (PSpos(1) + SSpos(1))/2 + offsetX, (PSpos(2) + SSpos(2))/2 + offsetY, '$\distLinePrimSource$');
t(5) = text(ax, (RecPos(1) + SSpos(1))/2 + offsetX, (RecPos(2) + SSpos(2))/2 + offsetY, '$\distLinePoint$');
t(6) = text(ax, posMidArc1(1), posMidArc1(2) + 0.01, '$\normPrimaryPropAngleSection$');
t(7) = text(ax, posMidArc2(1), posMidArc2(2) - 0.01, '$\normSecondPropAngleSection$');
t(8) = text(ax, SSpos(1) + longArrow*cosd(normalDir), SSpos(2) + longArrow*sind(normalDir), '$\surfaceNormal$');
t(9) = text(ax, 0.1, 0.8, '$\wfsPlane$');
t(10) = text(ax, 0.8, 0.1, '$\sectionTheo$');

% Adjust positions manually with the GUI

% % Then, reset the z position to 0
% for k = 1:numel(t)
%     t(k).Position(3) = 0;
% end

% % In the GUI, edit the text but without modify it just to "remind" matlab
% % that they are text boxes. I don't know why, but if you don't do it, it
% % doesn't work.

% % Transparency
% s.AlphaDataMapping = 'none';
% s.FaceAlpha = 'flat';
% s.AlphaData = 0.4*ones(2, length(xSect));
% faceNorm = permute(s.FaceNormals, [2 3 1]);
% % ax.View = [-20, 55];
% [camPos] = pointWiseExtend(ax.CameraPosition, faceNorm);
% visibleFaces = dot(faceNorm, camPos, 2) > 1;
% s.AlphaData(1:2, ~visibleFaces) = 0.4;
% s.AlphaData(1:2, visibleFaces) = 0;
% 
% currentFolder = pwd;
% cd(imagesPath); % Needed for inkscape to link svg files properly
% Plot2LaTeX( ax.Parent, 'pruebaTheoScheme'); % 'Kirchoff2_5DTheoScheme'
% cd(currentFolder)


%% Kirchhoff Integral

ax = axes(figure, 'NextPlot', 'Add', 'DataAspectRatio', [1 1 1]);

% Create volume of integration
[Vx, Vy, Vz] = ellipsoid(0, 0, 0, 2, 1, 1, 30);
col = [0.9490    0.9020    0.3765];
colMat = repmat(permute(col, [1 3 2]), [size(Vx), 1]);
s = surf(ax, Vx, Vy, Vz, colMat);
s.LineStyle = 'none';

% Lighting properties
s.FaceLighting = 'gouraud';
s.BackFaceLighting = 'lit';
% s.FaceNormals = -s.FaceNormals;
l = light(ax, 'Position', [-100 0 0], 'Style', 'local');

% Create primary source
PSpos = [-3 0 1];
ps = scatter3(ax, PSpos(1), PSpos(2), PSpos(3), 50, [0 0 0], 'filled');

% Create receiving point
RecPos = [0 0 0];
recP = scatter3(ax, RecPos(1), RecPos(2), RecPos(3), 50, [0 0 0], 'filled');

% Create normal arrow to S in section at a random point
theta = 180; height = 0.9;
[Vtheta, rho, Vheight] = cart2pol(Vx, Vy, Vz);
% Find the current vertex with the closest orientation to the ideally
% defined by theta and height
aux = [cos(Vtheta(:)), sin(Vtheta(:)), Vheight(:)];
dist = calcDistances(aux, [cosd(theta), sind(theta), height]);
[~, ind] = min(dist);
SSpos = [Vx(ind), Vy(ind), Vz(ind)]; % Secondary source position
[xInd, yInd] = ind2sub(size(Vx), ind);
normalDir = -squeeze(s.FaceNormals(xInd, yInd, :))';
longArrow = 0.5;
arrow(SSpos, SSpos + normalDir*longArrow)

ax.Visible = 'off';

% % In the GUI, edit the text but without modify it just to "remind" matlab
% % that they are text boxes. I don't know why, but if you don't do it, it
% % doesn't work.

% Transparency
s.AlphaDataMapping = 'none';
s.FaceAlpha = 'flat';
s.AlphaData = 0.2*ones(size(Vx));
% Spoints = [Vx(:), Vy(:), Vz(:)];
% faceNorm = Spoints./repmat(vecnorm(Spoints, 2, 2), 1, 3);
% % ax.View = [-20, 55];
% [x, y, z] = sph2cart(deg2rad(ax.View(1) - 90), deg2rad(ax.View(2)), 1);
% camDir = [x, y, z];
% [camDirMat] = pointWiseExtend(camDir, faceNorm);
% visibleFaces = dot(faceNorm, camDirMat, 2) > 0;
% s.AlphaData(~visibleFaces) = 0.4;
% s.AlphaData(visibleFaces) = 0.1;

text(ax, PSpos(1), PSpos(2), PSpos(3), '$\PosTheo[primarySource]$');
text(ax, SSpos(1) + longArrow*normalDir(1), SSpos(2) + longArrow*normalDir(2),...
    SSpos(3) + longArrow*normalDir(3), '$\surfaceNormal$');
text(ax, RecPos(1), RecPos(2), RecPos(3), '$\PosTheo$');
text(ax, SSpos(1), SSpos(2), SSpos(3), '$\PosTheo[surface]$');
text(ax, 2, 2, 2, '$V$')

% currentFolder = pwd;
% cd(imagesPath); % Needed for inkscape to link svg files properly
% Plot2LaTeX( ax.Parent, 'pruebaTheoScheme'); % 'KirchhoffTheoScheme'
% cd(currentFolder)

%% Rayleigh integrals scenario

ax = axes(figure, 'NextPlot', 'Add', 'DataAspectRatio', [1 1 1], ...
    'XLim', [-1 1], 'YLim', [-1 1], 'ZLim', [-0.5 1]);
ax.View = [-20, 10];

% Draw circle 
t = linspace(0, 2*pi);
r = 1;
x = r*cos(t);
y = r*sin(t);
figure(1)
S1 = patch(ax, x, y, 'yellow', 'FaceAlpha', 0.5);

% Draw hemisphere
theta = linspace(0, pi/2, 10)';
phi = linspace(0, 2*pi, 20);
x = cos(theta)*cos(phi);
y = cos(theta)*sin(phi);
z = pointWiseExtend(sin(theta), x);
S2 = surf(ax, x, y, z, 'FaceColor','yellow','FaceAlpha',0.5, 'LineStyle', 'none');

% Create primary source
PSpos = [0 0 -0.5];
ps = scatter3(ax, PSpos(1), PSpos(2), PSpos(3), 50, [0 0 0], 'filled');

% Create receiving point
RecPos = [0.5 0 0.5];
recP = scatter3(ax, RecPos(1), RecPos(2), RecPos(3), 50, [0 0 0], 'filled');

% Create secondary source point and arrow
SSpos = [0 0 0];

% Create normal arrow
normalDir = [0 0 1];
longArrow = 0.25;
arrow(SSpos, SSpos + normalDir*longArrow)
% 
% % Lighting properties
% S2.FaceLighting = 'gouraud';
% S2.BackFaceLighting = 'lit';
% % s.FaceNormals = -s.FaceNormals;
% l = light(ax, 'Position', [-100 0 0], 'Style', 'local');

% Transparency
S2.AlphaDataMapping = 'none';
S2.FaceAlpha = 'flat';
S2.AlphaData = 0.2*ones(size(x));

t = gobjects(0);
t = [t, text(ax, PSpos(1), PSpos(2), PSpos(3), '$\PosTheo[primarySource]$')];
t = [t, text(ax, SSpos(1) + longArrow*normalDir(1), SSpos(2) + longArrow*normalDir(2),...
    SSpos(3) + longArrow*normalDir(3), '$\surfaceNormal$')];
t = [t, text(ax, RecPos(1), RecPos(2), RecPos(3), '$\PosTheo$')];
t = [t, text(ax, SSpos(1), SSpos(2), SSpos(3), '$\PosTheo[surface]$')];
t = [t, text(ax, 2, 2, 0, '$\surfaceTheoRayleighSemisphere$')];
t = [t, text(ax, 0, 0, 0, '$\surfaceTheoRayleighPlane$')];

% for tcurr = t
%     outsideLimitsX = tcurr.Position(1) < ax.XLim(1) || tcurr.Position(1) > ax.XLim(2);
%     outsideLimitsY = tcurr.Position(2) < ax.YLim(1) || tcurr.Position(2) > ax.YLim(2);
%     outsideLimitsZ = tcurr.Position(3) < ax.ZLim(1) || tcurr.Position(3) > ax.ZLim(2);
%     if outsideLimitsX
%         tcurr.Position(1) = 0;
%     end
%     if outsideLimitsY
%         tcurr.Position(2) = 0;
%     end
%     if outsideLimitsZ
%         tcurr.Position(3) = 0;
%     end
% end

S2.Visible = 'off';
S1.Visible = 'off';
S2.Visible = 'on';
S1.Visible = 'on';

ax.Visible = 'off';

% currentFolder = pwd;
% cd(imagesPath); % Needed for inkscape to link svg files properly
% Plot2LaTeX( ax.Parent, 'pruebaTheoScheme'); % 'RayleighTheoScheme'
% cd(currentFolder)
