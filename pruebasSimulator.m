%% Integral de Kirchhoff para una fuente puntual y una superficie cúbica
% Source
r0 = [0 0 0];
k = 1;
A = 1;

% Representation Plane
XLim_plane = [0 80];
YLim_plane = [0 80];
XnumPoints_plane = 50;
YnumPoints_plane = 50;
r_plane = plane(XLim_plane, YLim_plane, XnumPoints_plane, YnumPoints_plane, [], [], []);

% Configure simulator with one source
objMono = simulator;
objMono.sourcePositions = r0;
objMono.sourceCoefficients = A;
objMono.k = k;

obj = simulator(figure);
obj.sourcePositions = r0;
obj.sourceCoefficients = A;
obj.sourceOrientations = [0 1 0 0];
obj.radPatFuns = {@(x) simulator.monopoleRadPat(x)};
obj.k = k;

obj.XLim = XLim_plane;
obj.YLim = YLim_plane;
obj.XnumPoints = XnumPoints_plane;
obj.YnumPoints = YnumPoints_plane;

% Animation parameters
t = 0:0.25:20;
f = 1;

% A) The real field
obj.simulate();
obj.animate(f, t);

U_orig = sphericalWave(k, A, r0, r_plane);
U_orig(abs(U_orig) > 2) = 2;
fig = singleFreqRepr( r_plane, U_orig, t, f, [-1, 1] );

% B) The new sources have the same complex coefficient that the signal that
% arrives to them from the original source

% Closed surface parameters
Lim = [20 40];
numPoints = 20;
XLim_surface = Lim;
YLim_surface = Lim;
ZLim_surface = [-10 10];
XnumPoints_surface = numPoints;
YnumPoints_surface = numPoints;
ZnumPoints_surface = 10;

[r_surface, dS_vec] = prism(XLim_surface, YLim_surface, ZLim_surface, XnumPoints_surface, YnumPoints_surface, ZnumPoints_surface);
numPointsSurface = size(r_surface, 1);
dS = modVec(dS_vec);
n = scaleVector(dS_vec, 1./dS);

% Huygens simulator
U_surface_simulator = objMono.calculate(r_surface);
obj.sourcePositions = r_surface;
obj.sourceOrientations = simulator.vec2rotVec(n);
obj.sourceCoefficients = dS.*U_surface_simulator;
obj.radPatFuns = repmat({@(~) 1}, numPointsSurface, 1);
obj.simulate();
obj.animate(1,t);

% Huygens hardcoded
U_surface = sphericalWave(k, A, r0, r_surface);
U_huygens = sphericalWave(k, dS.*U_surface, r_surface, r_plane);
fig = singleFreqRepr( r_plane, {U_orig, U_huygens}, t, f, [-0.1, 0.1] );

% C) We follow Kirchhoff's integral
% Source
r0 = [0 0 0];
k = 1;
A = 1;

% Closed surface parameters
Lim = [50 75];
numPoints = 50;
XLim_surface = Lim;
YLim_surface = Lim;
ZLim_surface = [50 150];
XnumPoints_surface = numPoints;
YnumPoints_surface = numPoints;
ZnumPoints_surface = 100;

[r_surface, dS_vec] = prism(XLim_surface, YLim_surface, ZLim_surface, XnumPoints_surface, YnumPoints_surface, ZnumPoints_surface);
numPointsSurface = size(r_surface, 1);
dS = modVec(dS_vec);
n = scaleVector(dS_vec, 1./dS);

% Measure point
p = [65, 60, 100];
U_p_real = sphericalWave(k, A, r0, p);

% Method 1
U_surface = sphericalWave(k, A, r0, r_surface);
% Kirchhoff's integral
r_r0 = r_surface - repmat(r0, numPointsSurface, 1);
grad_U = A*gradGreenFunction(r_r0, k);

[U1, U2] = kirchhoffIntegral(k, r_surface, dS_vec, U_surface, grad_U, p);

% Hardcoded
r = p;
U = U_surface;
area_dS = dS;

r_rs = repmat(r, numPointsSurface, 1) - r_surface;
r_rs_mod = modVec(r_rs);

G = exp(-1i*k*r_rs_mod)./r_rs_mod;
grad_G = scaleVector(-r_rs, exp(-1i*k*r_rs_mod)./r_rs_mod.^2.*(-1i*k - 1./r_rs_mod));

U_r = 1/(4*pi)*sum( area_dS.* (U.*dot(n, grad_G, 2) - G.*dot(n, grad_U, 2)) );
aux = area_dS.* (U.*dot(n, grad_G, 2) - G.*dot(n, grad_U, 2));
aux1 = U.*dot(n, grad_G, 2).*area_dS/(4*pi);
aux2 = - G.*dot(n, grad_U, 2).*area_dS/(4*pi);

% Method 2 (Simulator)
obj.setKirchhoffIntegralScenario(r0, 1, r_surface, -n, dS)
obj.simulate();
obj.animate(1,t);
obj.sourceCoefficients(1) = 0;

