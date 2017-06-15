%% Parámetros del escenario
% Definir los parámetros del escenario

% Global
c = 340; % m/s
freq = 440; % Frequency in Hz

% Noise source
noiseSourcesPosition = [0 0 0];
sourceCoef = 1;
numNoiseSources = size(noiseSourcesPosition, 1);

% Cancelation array
[ x, y, alfa ] = octogon(0.18, 8, 8, 24, 45);
loudspeakersPosition = [x, y, zeros(numel(x), 1)];
loudspeakersOrientation = [cosd(alfa), sind(alfa), zeros(numel(alfa), 1)];
numLouds = size(loudspeakersPosition, 1);

numSources = numLouds + numNoiseSources;

% Measure points
measurePos = [3 3 0];
numMeasurePoints = size(measurePos, 1);

%% Cálculo de envolventes que alimentan los altavoces
scenObj = scenario;
scenObj.setScenario(noiseSourcesPosition, loudspeakersPosition, loudspeakersOrientation);

delays = scenObj.delays;
attenuations = scenObj.attenuations;

loudspeakerCoef = repmat(sourceCoef.', [numLouds, 1]).*attenuations.*exp(1i*2*pi*freq*delays);
loudspeakerCoef = sum(loudspeakerCoef, 2);

%% Simulación
simulObj = simulator;

% Configure simulator object
simulObj.sourcePositions = [noiseSourcesPosition; loudspeakersPosition];
simulObj.sourceCoefficients = [sourceCoef; loudspeakerCoef];
simulObj.sourceOrientations = simulator.vec2rotVec([repmat([0 0 1], [numNoiseSources, 1]); loudspeakersOrientation]); % Mx4 matrix. Rotation vector: [angle of rotation, Xaxis, Yaxis, Zaxis]
simulObj.radPatFuns = repmat({@(x) simulator.monopoleRadPat(x)}, [numSources, 1]);
simulObj.k = 2*pi*freq/c;

% Simulate
U = simulObj.calculate(measurePos);

%% Reproducción y grabación
reprodObj = reproductorRecorder;

% Configure reproductorRecorder object

% Execute reproduction

% Get and analyze recorded signals
recorded = reprodObj.recorded;



%% Comparación de simulación y mediciones





    
            
                    