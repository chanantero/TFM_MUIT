%% Tests to see if the reproduction and recording in the laboratory works

%% A) Reproduction
sampleRate = 44100;
deviceWriter = audioDeviceWriter;
deviceWriter.SampleRate = sampleRate;
deviceWriter.Driver = 'ASIO';
deviceWriter.Device = 'Default';
deviceWriter.ChannelMappingSource = 'Property';

signProv = signalProvider;
filename = 'C:\Users\Rubén\Music\Salsa\Flor Pálida - Marc Anthony.mp3';
signProv.FileName = filename;
signProv.SamplesPerFrame = sampleRate*2; % The maximum an audioDeviceWriter allows

tic
while ~isDone(signProv)
    x = step(signProv);
    step(deviceWriter, x(:, [1 2]));
end
toc
disp('done')

release(signProv)
release(deviceWriter)

%% B) Reproduction and recording with ReproductorRecorder
% On 30/08/2018, the class ReproductorRecorder does not include yet the use
% of audioPlayerRecorder

repRecObj = reproductorRecorder;
repRecObj.setProps('mode', originType('custom'), 1);
repRecObj.setProps('enableProc', false);
repRecObj.setProps('Fs_player', sampleRate, 1);
repRecObj.setProps('Fs_recorder', sampleRate, 1);

% - Only noise source
aux = customSignal;
aux(:, wfsIndDest) = 0;
repRecObj.setProps('customSignal', aux, 1);
repRecObj.executeOrder('play');
recSignalNS = repRecObj.recorded{1}.';

% - Only WFS
aux = customSignal;
aux(:, nsIndDest) = 0;
repRecObj.setProps('customSignal', aux);
repRecObj.executeOrder('play');
recSignalWFS = repRecObj.recorded{1}.';

% - All
repRecObj.setProps('customSignal', customSignal);
repRecObj.executeOrder('play');
recSignal = repRecObj.recorded{1}.';