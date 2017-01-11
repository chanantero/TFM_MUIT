%% Full input loaded and full delay and attenuation known
% Variables needed for generating output signal.
% - input signal
% - sampling frequency
% - delay matrix
% - attenuation matrix.

fileName = 'Alan Walker - Faded.mp3';

fileInfo = audioinfo(fileName);
Fs = fileInfo.SampleRate;
input = audioread(fileName);
input = mean(input, 2); % Average of all input channels. Mono signal

% Apply delay
delayed = applyDelay( input, Fs, delays );

% Apply attenuation
output = delayed.*attenuations;

% Reproduce in loops.
frameSize = 1024*8;
output = input;
reproduceLoops( output, Fs, frameSize );


%% Reproduce with certain values of attenuation and delay

numChannels = 2; % Number of output channels
frameSize = 1024;

fileName = '02 - Dangerous Woman.mp3';

fileReader = dsp.AudioFileReader(fileName , 'SamplesPerFrame', frameSize);
fileInfo = audioinfo(fileName);

deviceWriter = audioDeviceWriter('SampleRate', fileInfo.SampleRate);
setup(deviceWriter, zeros(fileReader.SamplesPerFrame,fileInfo.NumChannels));

frameCount = 0;
while ~isDone(fileReader) && frameCount < 20
    frameCount = frameCount + 1;
    audioInput = step(fileReader);
    delays = ; % Vector with numChannels elements. Each one is the delay
    attenuations = ; % Ídem
    play(deviceWriter, audioInput);
end

%%
% Close the input file and release the device. The reproduction stops!!
release(fileReader);
release(deviceWriter);

%% With complete audio file import

% Input file
fileName = 'Alan Walker - Faded.mp3';

fileInfo = audioinfo(fileName);
audioInput = audioread(fileName, [1, fileInfo.SampleRate*10]); %[1, fileInfo.TotalSamples]);
audioInput = mean(audioInput, 2); % Average of all input channels. Mono signal
numSamples = numel(audioInput);

% Processing parameters
numChannels = 2; % Number of output channels

delays = zeros(numSamples, numChannels); % numSamples x numChannels array. 
% The (i,j) element corresponds to the delay of the signal at the moment 
% (i-1)/SampleFrequency in the j-th channel.

attenuation = ones(numSamples, numChannels); % numSamples x numChannels 
% array. The (i,j) element corresponds to the attenuation of the signal at
% the moment (i-1)/SampleFrequency in the j-th channel.

% Simulate a change of position in a circunference
w = pi/2; % rad/s
r = 2.5; % meters
c = 340; % meters/s
sepEars = 0.2; % Separation between ears
posRightEar = [sepEars/2, 0];
posLeftEar = [-sepEars/2, 0];

t = (0:numSamples-1)'/fileInfo.SampleRate;
posX = cos(w*t)*r;
posY = sin(w*t)*r;

distRightEar = [posX posY] - repmat(posRightEar, numel(t), 1);
distRightEar = sqrt(sum(distRightEar.^2, 2));
distLeftEar = [posX posY] - repmat(posLeftEar, numel(t), 1);
distLeftEar = sqrt(sum(distLeftEar.^2, 2));

delaysRightEar = distRightEar/c;
delaysLeftEar = distLeftEar/c;
delays = [delaysLeftEar, delaysRightEar];

attenuationRightEar = 1./distRightEar;
attenuationLeftEar = 1./distLeftEar;
attenuation = [attenuationLeftEar, attenuationRightEar];

% Processing
delays_sampleUnits = round(delays*fileInfo.SampleRate);
indices = repmat((1:numSamples)', 1, numChannels) - delays_sampleUnits;

audioOutput = repmat(audioInput, 1, numChannels);
audioOutput(indices > 0) = audioInput(indices(indices > 0));
audioOutput(indices <= 0) = 0;
audioOutput = audioOutput.*attenuation;

% Output streaming

% Send all the audioOutput to the buffer at once. Be careful,  the memory
% could be not big enough
% deviceWriter = audioDeviceWriter('SampleRate', fileInfo.SampleRate, 'SupportVariableSizeInput', true);
setup(deviceWriter, zeros(numSamples, numChannels));
% play(deviceWriter, audioOutput);
play(deviceWriter, zeros(10000, 2))

% In loops. More scalable, less memory requirements
frameSize = 1024*8;
deviceWriter = audioDeviceWriter('SampleRate', fileInfo.SampleRate);
setup(deviceWriter, zeros(frameSize, numChannels));

% For loop
numFrames = floor(numSamples/frameSize);
extraSamples = rem(numSamples,frameSize);
for k = 1:numFrames
    play(deviceWriter, audioOutput((k-1)*frameSize+1:k*frameSize, :));
end
if extraSamples>0
    play(deviceWriter, [audioOutput(numFrames*frameSize+1:end, :); zeros(frameSize - extraSamples, numChannels)]);
end

% While loop
k = 0;
finished = false;
while ~finished
    k = k+1;
    if k*frameSize >= numSamples
        play(deviceWriter, [audioOutput((k-1)*frameSize+1:end, :); zeros(k*frameSize - numSamples,numChannels)]);
        finished = true;
    else
        play(deviceWriter, audioOutput((k-1)*frameSize+1:k*frameSize, :));
    end  
end

% Release
release(deviceWriter);
