obj = simulator(figure);

% obj.simulate();

sourcePos = [0 0 0];
sourceCoeff = 1;

% Closed surface
XLim = [5 10]; YLim = [5 10]; ZLim = [0 0];
XnumPoints = 5; YnumPoints = 5; ZnumPoints = 2;
[surfPos, surfNormal] = prism(XLim, YLim, ZLim, XnumPoints, YnumPoints, ZnumPoints);
numPointsSurface = size(surfPos, 1);
surfArea = modVec(surfNormal);
n = scaleVector(surfNormal, 1./surfArea);

obj.setKirchhoffIntegralScenario(sourcePos, sourceCoeff, surfPos, n, surfArea);
obj.k = 5;

obj.XLim = [0, 11];
obj.YLim = [0, 11];
obj.XnumPoints = 50;
obj.YnumPoints = 50;
obj.z = 0;

obj.simulate();
obj.ax.CLim = [0 0.5];
% scat = findobj(obj.ax, 'Type', 'scatter');
% delete(scat)

%%
obj = simulator(figure);
obj.k = 5;

sourcePos = [0 0 0];
sourceCoeff = 1;
surfacePointsPos = square([5 10], [5 10], 20, 20, [], [], []);


obj.setLinearArrayScenario(sourcePos, sourceCoeff, surfacePointsPos)

obj.XLim = [0, 11];
obj.YLim = [0, 11];
obj.XnumPoints = 100;
obj.YnumPoints = 100;
obj.z = 0;

obj.simulate();
obj.sourceCoefficients(1) = 0;
obj.sourceCoefficients(1) = 1;

obj.ax.CLim = [0 1];


%% Only one source
obj = simulator(figure);
obj.ax.CLim = [0 1];

sourcePos = [0 0 0];
sourceCoeff = 1;

obj.setSources(sourcePos, 'coefficient', sourceCoeff);

obj.XLim = [0, 5];
obj.YLim = [0, 5];
obj.XnumPoints = 100;
obj.YnumPoints = 100;
obj.z = 0;

obj.k = 5;

obj.simulate();
f = 1;
t = 0:0.1:5;
obj.animate(f, t);


%% One line
obj = simulator(figure);
obj.k = 2;

sourcePos = [0 0 0];
sourceCoeff = 1;
x = (-50:1:50)';
y = 5*ones(numel(x), 1);
z = zeros(numel(x), 1);
linePointsPos = [x, y, z];

obj.setLinearArrayScenario(sourcePos, sourceCoeff, linePointsPos)

obj.XLim = [-60, 60];
obj.YLim = [0, 11];
obj.XnumPoints = 300;
obj.YnumPoints = 300;
obj.z = 0;

obj.simulate();
f = 1;
t = 0:0.1:5;
obj.animate(f, t);

obj.sourceCoefficients(1) = 0;
obj.sourceCoefficients(1) = 1;
obj.sourceCoefficients(2:end) = 0;

obj.ax.CLim = [0 0.1];

