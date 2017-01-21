classdef WFSTool2 < handle
    
    properties
        fig
        player
    end
    
    methods
        
        function obj = WFSTool2()
            obj.fig = openfig('WFSTool.fig');
            
            Fs = 44100;
            numOutputChannels = 2;

            % Reading object
            obj.fileReader = dsp.AudioFileReader();
            
            % Processing object
            obj.processor = processSignal('Fs', Fs, 'variable', true, 'numChannels', numOutputChannels, 'delayType', 'forward');
            
            % Writing object
            frameSizePlaying = Fs;
            obj.player = audioPlayer('Fs', Fs, 'numChannels', numOutputChannels, 'frameSize', frameSizePlaying);
            
            % Modify original callbacks
            butTest = findobj(obj.fig.Children, 'Tag', 'but_test');
            butTest.Callback = @(hObject, eventdata) obj.callback(); %@(hObject, eventdata) but_test_Callback(obj.playObj, Fs);
            
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
            obj.fileReader.Filename = fileName;
                               
            % Reproduce
            frameCount = 0;
            while ~isDone(obj.fileReader) %&& frameCount < 20
                frameCount = frameCount + 1;
                delay = zeros(numOutputChannels, 2);
                step(player, delay);
                pause(0.001)
            end
            release(obj.fileReader);
            release(obj.processor);
            release(obj.player);
            fprintf('Finished\n')
        end
       
        function executeOrder(obj, order)
            switch state
                case 'playing'
                    switch order
                        case 'stop'
                            obj.s
                        case 'play'
                        case 'resume'
                        case 'pause'
                            
                    end
                    
                case 'stopped'
                    
                case 'paused'
                    
            end
        end
        
    end
    
end

