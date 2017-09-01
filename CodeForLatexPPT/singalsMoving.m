% Animation of signals progressing in time
pulseCoefMat(:, :, 1) = ...
    [1 0;
    0 1;
    0 0;
    -1 0];

pulseCoefMat(:, :, 2) = ...
    [0 0;
    0.5 0;
    1i 1;
    0 1];

pulseLimits = 1*...
    [1 2;
    3 4;
    5 7;
    9 10];
freqs = [2; 4];
sampleRate = 10000;

% Channel impulse response
h1 = 1;
h2 = 0.5*1i;

% Graphics
fig1 = figure;
fig2 = figure;
fig3 = figure;
ax1 = axes(fig1);
ax2 = axes(fig2);
ax3 = axes(fig3);
YLim = [-1 1];

% Animation parameters
windowSpan = 4;
framePerSecond = 30; % 1/Time-step
t_ini = -windowSpan:1/framePerSecond:pulseLimits(end, 2)+windowSpan;
numFrames = numel(t_ini);

% Saving parameteres
path = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\';
fileName1 = 'movingWaveTrans1.gif';
fileName2 = 'movingWaveTrans2.gif';
fileName3 = 'movingWaveRec.gif';

% Initialize
t0 = t_ini(1);
t_end = t0 + windowSpan;
signal = pulseCoefMat2signal(pulseCoefMat, pulseLimits, freqs, sampleRate, t0, t_end, 'type_marker', 'time', 'type_pulseLimits', 'time');
numSamples = size(signal, 1);
t = (0:numSamples-1)/sampleRate;
l1 = plot(ax1, t, signal(:, 1));
l1.LineWidth = 2;
l1.Color = [0 176 240]/255;
ax1.YLim = YLim;
ax1.XTick = [];
ax1.YTick = [];
F = getframe(ax1);
[GIFframe1, map1] = rgb2ind(F.cdata, 64);
imwrite(GIFframe1, map1, [path, fileName1],...
    'LoopCount', 100, 'DelayTime', 1/framePerSecond);


l2 = plot(ax2, t, signal(:, 2));
l2.LineWidth = 2;
l2.Color = [255, 192, 0]/255;
ax2.YLim = YLim;
ax2.XTick = [];
ax2.YTick = [];
F = getframe(ax2);
[GIFframe2, map2] = rgb2ind(F.cdata, 64);
imwrite(GIFframe2, map2, [path, fileName2],...
    'LoopCount', 100, 'DelayTime', 1/framePerSecond);


l3 = plot(ax3, t, signal);
arrayfun(@(l) eval('l.LineWidth = 2;'), l3);
l3(1).Color = [0 176 240]/255;
l3(2).Color = [255, 192, 0]/255;
l3(1).YData = 0.5*ones(size(l3(1).XData));
ax3.YLim = YLim;
ax3.XTick = [];
ax3.YTick = [];
F = getframe(ax3);
[~, map3] = rgb2ind(F.cdata, 64);
l3(1).YData = signal(:,1);
F = getframe(ax3);
GIFframe3 = rgb2ind(F.cdata, map3);
imwrite(GIFframe3, map3, [path, fileName3],...
    'LoopCount', 100, 'DelayTime', 1/framePerSecond);

% Process
for frame = 1:numFrames
t0 = t_ini(frame);
t_end = t0 + windowSpan;
signal = pulseCoefMat2signal(pulseCoefMat, pulseLimits, freqs, sampleRate, t0, t_end, 'type_marker', 'time', 'type_pulseLimits', 'time');
numSamples = size(signal, 1);
t = (0:numSamples-1)/sampleRate;

l1.XData = t;
l1.YData = signal(:,1);
F = getframe(ax1);
GIFframe = rgb2ind(F.cdata, map1);
imwrite(GIFframe, map1, [path, fileName1], 'WriteMode', 'append',...
    'DelayTime', 1/framePerSecond);

l2.XData = t;
l2.YData = signal(:,2);
F = getframe(ax2);
GIFframe = rgb2ind(F.cdata, map2);
imwrite(GIFframe, map2, [path, fileName2], 'WriteMode', 'append',...
    'DelayTime', 1/framePerSecond);

t0 = t_ini(frame) - windowSpan;
t_end = t_ini(frame);
signal = pulseCoefMat2signal(pulseCoefMat, pulseLimits, freqs, sampleRate, t0, t_end, 'type_marker', 'time', 'type_pulseLimits', 'time');
numSamples = size(signal, 1);
t = (0:numSamples-1)/sampleRate;
for indLine = 1:numel(l3)
    l3(indLine).XData = t;
    l3(indLine).YData = signal(:, indLine);
end
F = getframe(ax3);
GIFframe = rgb2ind(F.cdata, map3);
imwrite(GIFframe, map3, [path, fileName3], 'WriteMode', 'append',...
    'DelayTime', 1/framePerSecond);


end
