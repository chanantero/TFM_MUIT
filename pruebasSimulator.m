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
r_r0 = r_surface - repmat(r0, numPointsSurface, 1);
grad_U = sphericalWave_grad( k, A, r0, r_surface );
[U, A, B] = kirchhoffIntegral(k, r_surface, dS_vec, U_surface, grad_U, p);

% Method 2 (Simulator)
obj = simulator;
obj.setKirchhoffIntegralScenario(r0, 1, r_surface, n, dS)
U = obj.calculate(p);

% Hardcoded
measurePointPositions = p;

r_measure = p;
dS = dS_vec;
Usurf = U_surface;
grad_Usurf = grad_U;

numMeasurePoints = size(r_measure, 1);
numPointsSurface = size(r_surface, 1);
area_dS = modVec(dS);
n = scaleVector(dS, 1./area_dS);

U = zeros(numMeasurePoints, 1);
for l = 1:numMeasurePoints
    r = r_measure(l, :);
    rel_r = repmat(r, numPointsSurface, 1) - r_surface;
    
    G = GreenFunction(rel_r, k);
    grad_G = gradGreenFunction(-rel_r, k); % The minus sign is because we must find the derivate of (x - x0)
    A = G .* sum(n.*grad_Usurf, 2);
    B = Usurf .* sum(n.*grad_G, 2);
    U(l) = 1/(4*pi)*sum( area_dS.* (A - B) );
end

s = numPointsSurface + 1 + 1;
sourceCoef*4*pi/area_dS(s-1 -numPointsSurface)
Usurf(1)

G(s-numMeasurePoints)*radPatCoef
sum(n(s-25001,:).*grad_G(s-25001, :), 2)
U_aux(s)

plot(imag(A/(4*pi).*area_dS./U_aux((2:1+numPointsSurface)).'))
plot(real(B/(4*pi).*area_dS./U_aux((2:1+numPointsSurface)+numPointsSurface).'))


%% Sencilla

% Source
r0 = [0 0 0];
k = 1;
A = 1;

% Closed surface parameters
Lim = [25 50];
numPoints = 100;
XLim_surface = Lim;
YLim_surface = Lim;
ZLim_surface = [-10 10];
XnumPoints_surface = numPoints;
YnumPoints_surface = numPoints;
ZnumPoints_surface = 50;

[r_surface, dS_vec] = prism(XLim_surface, YLim_surface, ZLim_surface, XnumPoints_surface, YnumPoints_surface, ZnumPoints_surface);
numPointsSurface = size(r_surface, 1);
dS = modVec(dS_vec);
n = scaleVector(dS_vec, 1./dS);

% Real measure
p = [35, 35, 0];
U_p_real = sphericalWave(k, A, r0, p);

% Representation Plane
XLim_plane = [20 55];
YLim_plane = [20 55];
XnumPoints_plane = 20;
YnumPoints_plane = 20;

obj = simulator(figure);
obj.XLim = XLim_plane;
obj.YLim = YLim_plane;
obj.XnumPoints = XnumPoints_plane;
obj.YnumPoints = YnumPoints_plane;
obj.setKirchhoffIntegralScenario(r0, 1, r_surface, n, dS)

U = obj.calculate(p);

obj.simulate()
t = 0:0.25:20;
f = 1;
obj.animate(f, t)

obj.sourceCoefficients(1) = 1;