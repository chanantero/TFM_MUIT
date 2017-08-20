%% Change delay continuously
% Thought for single frequency tones
freq = 440;
T2pi = 2; % Time to do a phase shift of 2*pi radians
t = (0:0.01:60)';

N = numel(t);
delays = repmat(obj.scenarioObj.delays(:, 1)', N, 1);
attenuations = repmat(obj.scenarioObj.attenuations(:, 1)', N, 1);

channel = 48; % Independent channel of the noise source
delays(:, channel) = t./(freq*T2pi);
attenuations(:, channel) = 1;

obj.player.setProps('modeProc', timeInteractionTypes('realTime'), [1 1]);
obj.player.setProps('pred_t', t, [1 1]);
obj.player.setProps('pred_delay', delays, [1 1]);
obj.player.setProps('pred_atten', attenuations, [1 1]);


%%
% obj = WFSTool2_noiseChannel;

R = 5; % Radius
centerX = 1.7382; centerY = 3.1782;
t = (0:0.01:60)';
f = 0.1;
x = centerX + R*cos(2*pi*f*t);
y = centerY + R*sin(2*pi*f*t);
positions = [x, y, zeros(numel(t), 1)];

% From the positions, calculate the delays and attenuations using the
% scenario object
N = size(positions, 1);
numLoudspeakers = obj.scenarioObj.numLoudspeakers;
delays = zeros(N, numLoudspeakers);
attenuations = zeros(N, numLoudspeakers);
for k = 1:size(positions, 1)
    obj.scenarioObj.setSourcePosition(positions(k, :), 1);
    delays(k, :) = obj.scenarioObj.delays(:, 1)';
    attenuations(k, :) = obj.scenarioObj.attenuations(:, 1)';
end

% Assign the delays and attenuations to the reproductor object and change
% the mode of time interaction to predefined
% obj.player.setProps('modeProc', timeInteractionTypes('predefined'), [1 1]);

obj.player.setProps('modeProc', timeInteractionTypes('realTime'), [1 1]);
obj.player.setProps('pred_t', t, [1 1]);
obj.player.setProps('pred_delay', delays, [1 1]);
obj.player.setProps('pred_atten', attenuations, [1 1]);

% Now you can reproduce


%%
t = (0:0.1:60)';
f = 0.2;
x = 0.3*cos(2*pi*f*t);
y = 0.3*sin(2*pi*f*t);
positions = [x, y, zeros(numel(t), 1)];

% From the positions, calculate the delays and attenuations using the
% scenario object
loudspeakersPosition = [-0.1 0 0; 0.1 0 0];
loudspeakersOrientation = [1 0 0; -1 0 0];
roomPosition = [-1 -1 2 2];

fig = figure;
scen = scenario(fig);
scen.setScenario(positions(1, :), loudspeakersPosition, loudspeakersOrientation, roomPosition)

N = size(positions, 1);
numLoudspeakers = size(loudspeakersPosition, 1);
delays = zeros(N, numLoudspeakers);
attenuations = zeros(N, numLoudspeakers);
for k = 1:size(positions, 1)
    scen.setSourcePosition(positions(k, :), 1);
    delays(k, :) = scen.delays';
    attenuations(k, :) = scen.attenuations';
end

% Create object
obj = WFSTool2;

% Assign the delays and attenuations to the reproductor object and change
% the mode of time interaction to predefined
obj.player.setProps('modeProc', timeInteractionTypes('predefined'), [1 1]);
obj.player.setProps('pred_t', t, [1 1]);
obj.player.setProps('pred_delay', delays, [1 1]);
obj.player.setProps('pred_atten', attenuations, [1 1]);

% Now you can reproduce


