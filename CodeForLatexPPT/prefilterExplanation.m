duration = 10;
sampleRate = 100;
numSamp = floor(duration*sampleRate);
t = (0:numSamp - 1)/sampleRate;

durPulse = 0.5;
tStart = 1;
f = 1/durPulse;
tEnd = tStart + durPulse;
ind = t >= tStart & t < tEnd;
x = zeros(numSamp, 1);
x(ind) = cos(2*pi*f*(t(ind) - tStart) - pi/2);

% ax = axes(figure);
% plot(ax, t, x)

tStart_h = -1; tEnd_h = 1;
numSamp_h = floor(tEnd_h*sampleRate) - floor(tStart_h*sampleRate) + 1;
t_h = (0:numSamp_h - 1)/sampleRate + tStart_h;

h = zeros(numSamp_h, 1);
ind = t_h >= -1/2 & t_h < 0;
h(ind) = -(1 + 2*t_h(ind));
ind = t_h >= 0 & t_h <= 1/2;
h(ind) = (1 - 2*t_h(ind));

% ax = axes(figure);
% plot(ax, t_h, h)

y_delayed = filter(h, 1, x)/sampleRate;
t_y_delayed = t - (-tStart_h - 1/2);

% ax = axes(figure);
% plot(ax, t, y_delayed);

% y = y_delayed;
% t_y = t + tStart_h;
% plot(ax, t_y, y)

x2 = x/2;
tau = 1.5;
t_x2 = t + tau;

axFilt = axes(figure);
plot(axFilt, t_h, h)
axFilt.XTick = [-1/2, 0, 1/2];
axFilt.XTickLabels = {'-T/2', '0', 'T/2'};
axFilt.YTick = [];
% axFilt.XLabel.String = 'Tiempo (s)';
% axFilt.YLabel.String = 'h(t)';

ax = axes(figure, 'NextPlot', 'Add');
plot(ax, t, x)
plot(ax, t_y_delayed, y_delayed, 'Color', [255, 192, 0]/255)
% plot(ax, t_y, y)
plot(ax, t_x2, x2, 'Color', [237, 125, 49]/255)
plot(ax, t_y_delayed + (tau - 1/2), y_delayed, 'Color', [0, 176, 80]/255)
plot(ax, [0, 10], [0, 0], 'k')

ax.XLim = [0.5, 3.5];
ax.XTick = [];
ax.YTick = [];
ax.Box = 'on';