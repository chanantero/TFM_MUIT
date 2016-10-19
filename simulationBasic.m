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
rSources = [0 0 0; 1 1 0; 0 1 0; 2 0 0];
XnumPoints = 100; YnumPoints = 100;
xVec = linspace(-1, 4, XnumPoints)'; yVec = linspace(-1, 2, YnumPoints)';
rImage = [kron(xVec, ones(YnumPoints, 1)), repmat(yVec, XnumPoints, 1), zeros(XnumPoints*YnumPoints, 1)];

U = sum(sphericalWave(k, ones(size(rSources, 1), 1), rSources, rImage), 2);

U = reshape(U, YnumPoints, XnumPoints);

U(real(U) > 10) = 10;

image(xVec, yVec, abs(U), 'CDataMapping', 'scaled')

%% Integral de Kirchhoff para una fuente puntual y una superficie cúbica
% Source
r0 = [0 0 0];
k = 1;
A = 1;

% Closed surface
Lim = [50 75];
numPoints = 100;
[r_surface, dS_vec] = prism(Lim, Lim, [-0.5, 0.5], numPoints, numPoints, 10);
numPointsSurface = size(r_surface, 1);
dS = modVec(dS_vec);
n = scaleVector(dS_vec, 1./dS);
U_surface = sphericalWave(k, A, r0, r_surface);

% Representation parameters
r_plane = plane([0, 80], [0, 80], 200, 200, [], [], []);
t = 0:0.25:20;
f = 1;

% Each element of the surface is a new isotropic source. What complex coefficient
% defines the signal of each new source?

% A) The real field
U_orig = sphericalWave(k, A, r0, r_plane);
U_orig(abs(U_orig) > 2) = 2;

fig = singleFreqRepr( r_plane, U_orig, t, f, [-1, 1] );


fig = singleFreqRepr( r_plane, {U_orig, U_huygens}, t, f, [-0.1, 0.1] );

scatter3(r_plane(:,1), r_plane(:,2), r_plane(:,3), [], abs(U_huygens./U_orig), '.');
colorbar

% B) The new sources have the same complex coefficient that the signal that
% arrives to them from the original source
U_huygens = sphericalWave(k, dS.*U_surface, r_surface, r_plane);

fig = singleFreqRepr( r_plane, U_huygens, t, f, [-1, 1] );

% C) We follow Kirchhoff's integral
p = [60, 60, 100];
p2 = [70, 70, 0];

r_p = r_surface - repmat(p, numPointsSurface, 1);
r_r0 = r_surface - repmat(r0, numPointsSurface, 1);
r_r0_mod = modVec(r_r0);

U_surface_virtual = 1/(4*pi)*dS.*U_surface.*...
    (- 1i*k*dot(scaleVector(r_p, 1./modVec(r_p)) - scaleVector(r_r0, 1./modVec(r_r0)), n, 2) ...
    - dot(scaleVector(r_p, modVec(r_p).^(-2)), n, 2) ...
    + dot(scaleVector(r_r0, modVec(r_r0).^(-2)), n, 2));

U_p = sphericalWave(k, U_surface_virtual, r_surface, p);
U_p2 = sphericalWave(k, U_surface_virtual, r_surface, p2);

U_p_real = sphericalWave(k, A, r0, p);
U_p2_real = sphericalWave(k, A, r0, p2);

U_p/U_p_real
U_p2/U_p2_real

% U = sphericalWave(k, U_surface_virtual, r_surface, r_plane);

% fig = singleFreqRepr( r_plane, U, t, f, [-0.01, 0.01] );

%% Probando la integral de Rayleigh
% Source
r0 = [0 0 0];
k = 1;
A = 1;

% Closed volume
Lim = [50 75];
numPoints = 50;
[r_surface, dS_vec] = prism(Lim, Lim, [50 150], numPoints, numPoints, 100);
numPointsSurface = size(r_surface, 1);
dS = modVec(dS_vec);
n = scaleVector(dS_vec, 1./dS);
U_surface = sphericalWave(k, A, r0, r_surface);


p = [60, 60, 100];

% Kirchhoff's integral
r_p = r_surface - repmat(p, numPointsSurface, 1);
r_r0 = r_surface - repmat(r0, numPointsSurface, 1);
r_r0_mod = modVec(r_r0);

grad_U = A*scaleVector(r_r0, exp(-1i*k*r_r0_mod)./r_r0_mod.^2.*(-1i*k - 1./r_r0_mod));

U_p_kirchhoff = kirchhoffIntegral(k, r_surface, dS_vec, U_surface, grad_U, p);

U_p_real = sphericalWave(k, A, r0, p);

% Rayleigh integral
[r_infPlane, n] = plane([-1000, 1000], [-100, 100], 10000, 10000, [], [], [0, 0, 50]);
numPointsInfPlane = size(r_infPlane, 1);
dS_infPlane = repmat(n, numPointsInfPlane, 1);
rPlane_r0 = r_infPlane - repmat(r0, numPointsInfPlane, 1);
rPlane_r0_mod = modVec(rPlane_r0);
grad_U = A*scaleVector(rPlane_r0, exp(-1i*k*rPlane_r0_mod)./rPlane_r0_mod.^2.*(-1i*k - 1./rPlane_r0_mod));

U_p_rayleigh = rayleighI(k, r_infPlane, dS_infPlane, grad_U, p);
U_p_rayleighVolume = rayleighI(k, r_surface, dS_vec, grad_U, p);
U = sphericalWave(k, A, r0, r_infPlane);
U_p_rayleighII = rayleighII(k, r_infPlane, dS_infPlane, U, p);

U_p_rayleigh/U_p_real
U_p_rayleighII/U_p_real
U_p_rayleighII/U_p_rayleigh
U_p_rayleigh/U_p_kirchhoff

