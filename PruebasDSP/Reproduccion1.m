%% Read from File and Write to Audio Device
% Read an MP3 audio file and play it through your default audio output
% device.

%%
% Create a |dsp.AudioFileReader| System object(TM) with default settings.
% Use the |audioinfo| function to return a structure containing information
% about the audio file.
fileReader = dsp.AudioFileReader('02 - Dangerous Woman.mp3',...
    'SamplesPerFrame', 1024*64);
fileInfo = audioinfo('02 - Dangerous Woman.mp3');

%%
% Create an |audioDeviceWriter| System object and specify the sample rate.
% Call |setup| to reduce the computational load of initialization in an
% audio stream loop.
deviceWriter = audioDeviceWriter(...
    'SampleRate',fileInfo.SampleRate);
setup(deviceWriter,...
    zeros(fileReader.SamplesPerFrame,fileInfo.NumChannels));

%%
% In an audio stream loop, read an audio signal frame from the file, and
% write the frame to your device.
index = 0;
while ~isDone(fileReader) && index < 20
    index = index + 1;
    audioData = step(fileReader);
%     if rem(index, 2)
%         audioData(:,1) = 0;
%     else
%         audioData(:,2) = 0;
%     end
    play(deviceWriter, audioData);
end

%%
% Close the input file and release the device. The reproduction stops!!
release(fileReader);
release(deviceWriter);

% fileReader = dsp.AudioFileReader('Alan Walker - Faded.mp3',...
%     'SamplesPerFrame', 1024);
% fileInfo = audioinfo('Alan Walker - Faded.mp3');
% 
% deviceWriter = audioDeviceWriter(...
%     'SampleRate',fileInfo.SampleRate);
% setup(deviceWriter,...
%     zeros(fileReader.SamplesPerFrame,fileInfo.NumChannels));
% 
% index = 0;
% while ~isDone(fileReader)
%     index = index + 1;
%     audioData = step(fileReader);
%     play(deviceWriter, audioData);
% end
