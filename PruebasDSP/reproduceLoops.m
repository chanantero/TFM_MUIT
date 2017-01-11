function reproduceLoops( signal, Fs, frameSize )
% Input arguments:
% - signal. Signal
% - Fs. Sampling frequency
% - frameSize. Size of the fragments of signal that will be sent to the
% audio driver

[numSamples, numChannels] = size(signal);

deviceWriter = audioDeviceWriter('SampleRate', Fs);
setup(deviceWriter, zeros(frameSize, numChannels));

k = 0;
finished = false;
while ~finished
    k = k+1;
    if k*frameSize >= numSamples
        sampleIndices = (k-1)*frameSize+1:numSamples;
        frame = [signal(sampleIndices); zeros(k*frameSize - numSamples, numChannels)]; % Fill with 0s
        finished = true;
    else
        sampleIndices = (k-1)*frameSize+1:k*frameSize;
        frame = signal(sampleIndices, :);
    end  
    play(deviceWriter, frame);
end

% Release
release(deviceWriter);

end

