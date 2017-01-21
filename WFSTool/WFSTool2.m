classdef WFSTool2 < handle
    
    properties
        fig
        fileReader
        processor
        playObj
    end
    
    methods
        
        function obj = WFSTool2()
            obj.fig = openfig('WFSTool.fig');
            
            Fs = 44100;
            numOutputChannels = 2;

            % Reading object
            frameSizeReading = Fs;
            obj.fileReader = dsp.AudioFileReader(fileName, 'SamplesPerFrame', frameSizeReading);
            
            % Processing object
            obj.procObj = processSignal('Fs', Fs, 'variable', true, 'numChannels', numOutputChannels, 'delayType', 'forward');
            
            % Writing object
            frameSizePlaying = Fs;
            obj.playObj = audioPlayer('Fs', Fs, 'numChannels', numOutputChannels, 'frameSize', frameSizePlaying);
            
            % Modify original callbacks
            butTest = findobj(obj.fig.Children, 'Tag', 'but_test');
            butTest.Callback = @(hObject, eventdata) but_test_Callback(obj.playObj, Fs);
            
            butImportWav = findobj(obj.fig.Children, 'Tag', 'file');
            table = findobj(obj.fig.Children, 'Tag', 'table');
            butImportWav.Callback = @(hObject, eventdata) open_Callback(table);
        end
        
        function playMusic(obj)
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
        end
        
    end
    
end

