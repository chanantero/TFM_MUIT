% Audio file information
fileName = '02 - Dangerous Woman.mp3';
fileInfo = audioinfo(fileName);

Fs = fileInfo.SampleRate;
numOutputChannels = 2;

% Reading object
frameSizeReading = Fs;
fileReader = dsp.AudioFileReader(fileName , 'SamplesPerFrame', frameSizeReading);

% Processing object
procObj = processSignal('Fs', Fs, 'variable', true, 'numChannels', numOutputChannels, 'delayType', 'forward');

% Writing object
frameSizePlaying = Fs;
playObj = audioPlaying('Fs', Fs, 'numChannels', numOutputChannels, 'frameSize', frameSizePlaying);

% Reproduce
frameCount = 0;
while ~isDone(fileReader) && frameCount < 20
    frameCount = frameCount + 1;
    audioInput = step(fileReader);
    audioInput = mean(audioInput, 2); % From Stereo to Mono
    delays = zeros(frameSizeReading, numOutputChannels);
    delays(:, 1) = 0.5;
    audioOutput = step(procObj, audioInput, delays);
    step(playObj, audioOutput);
end
fprintf('Finished\n')

% Finalize
release(fileReader)
release(procObj)
release(playObj)
