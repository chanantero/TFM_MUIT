% Spherical wave propagating
rSource = [0 0 0];
k = 1;

% We take a surface that encloses a volume
Lim = [1 2];
numPoints = 1000;
dS = (diff(Lim)/(numPoints-1))^2;
rSurface = prism(Lim, Lim, Lim, numPoints, numPoints, numPoints);

% Choose a point inside the volume, and compare the real value with the sum
% of inifinite (a lot) waves that come from diferential (very small)
% sources situated at the surface
rObj = [1.5, 1.5, 1.5]; % Point where we are going to measure and compare
U_ObjReal = sphericalWave(k, rSource, rObj); 

U_surface = sphericalWave(k, rSource, rSurface);
U_ObjVirtual = sum(conj(sphericalWave(k, rObj, rSurface)).*(U_surface*dS));

%% Convergencia según resolución del cubo
numPoints = round(logspace(log10(10), log10(1500), 30));
U_ObjVirtual = zeros(numel(numPoints), 1);
for k = 1:numel(numPoints)
    fprintf('k = %d\n', k)
    dS = (diff(Lim)/(numPoints(k)-1))^2;
    rSurface = prism(Lim, Lim, Lim, numPoints(k), numPoints(k), numPoints(k));
    U_surface = sphericalWave(k, rSource, rSurface);
    U_ObjVirtual(k) = sum(conj(sphericalWave(k, rObj, rSurface)).*(U_surface*dS));
end

ax = axes(figure, 'NextPlot', 'Add');
plot(ax, numPoints, abs(U_ObjVirtual))
% plot(ax, numPoints, angle(U_ObjVirtual))
ax.YScale = 'log';

%% Mapa
k = 2;
rSources = [0 0 0; 1 1 0; 0 1 0; 1 0 0];
XnumPoints = 100; YnumPoints = 100;
xVec = linspace(-1, 4, XnumPoints)'; yVec = linspace(-1, 2, YnumPoints)';
rImage = [kron(xVec, ones(YnumPoints, 1)), repmat(yVec, XnumPoints, 1), zeros(XnumPoints*YnumPoints, 1)];

U = sum(sphericalWave(k, rSources, rImage), 2);

U = reshape(U, YnumPoints, XnumPoints);

U(real(U) > 10) = 10;

image(xVec, yVec, abs(U), 'CDataMapping', 'scaled')