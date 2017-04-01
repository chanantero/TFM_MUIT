clear
obj = WFSTool2;
obj.scenarioObj.setForcedEnabledLoudspeakers([true, true]);

t = (0:0.1:60)';
f = 0.2;
x = cos(2*pi*f*t);
y = sin(2*pi*f*t);
positions = [x, y, zeros(numel(t), 1)];


obj.noRealTimeReproduction(t, positions);


%%
clear
obj1 = WFSTool2;
obj2 = WFSTool2;

obj1.player.setFrameSize(2*44100);
obj2.player.setFrameSize(2*44100);

t = (0:0.1:60)';
f = 0.2;
x = cos(2*pi*f*t);
y = sin(2*pi*f*t);
positions = [x, y, zeros(numel(t), 1)];

obj1.noRealTimeReproduction(t, positions);

