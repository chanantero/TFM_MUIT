classdef reproductor < matlab.System
     
    properties(Nontunable)
        frameSizeReading
        frameSizeWriting
        Fs
        audioFileName
    end
    
    properties(DiscreteState)

    end

    properties(Access = private)
        fileReader
        processor
        player
    end
    
    methods(Access = protected)
        
        function setupImpl(obj)
            
            
        end

        function stepImpl(obj, delay)
            
            % Read form file
            audioInput = step(obj.fileReader);
            audioInput = mean(audioInput, 2); % From Stereo to Mono
            
            % Process
            delays = repmat(delay', obj.frameSizeReading, 1);            
            audioOutput = step(obj.processor, audioInput, delays);
            
            % Write to audio device buffer
            step(obj.player, audioOutput);
            
        end

        function resetImpl(obj)
            obj.stop();
            obj.audioFileName = '';
        end

        function releaseImpl(obj)
            
        end
        
        
    end
    
    methods
        
        function obj = reproductor()
            obj.setDefaultProperties();

            % Reading object
            obj.fileReader = dsp.AudioFileReader();
            
            % Processing object
            obj.processor = processSignal('Fs', Fs, 'variable', true, 'numChannels', numOutputChannels, 'delayType', 'forward');
            
            % Writing object
            frameSizePlaying = Fs;
            obj.player = audioPlayer('Fs', Fs, 'numChannels', numOutputChannels, 'frameSize', frameSizePlaying);
            
        end
                
        function executeOrder(obj, order)
            switch state
                case 'playing'
                    switch order
                        case 'stop'
                            release(obj);
                        case 'play'
                        case 'resume'
                        case 'pause'
                            
                    end
                    
                case 'stopped'
                    
                case 'paused'
                    
            end
        end
        
        function stop(obj)
            release(obj.fileReader);
            release(obj.processor);
            release(obj.player);
        end
        
        function setDefaultProperties(obj)
            obj.frameSizeReading = 1024;
            obj.frameSizeWriting = 1024;
            obj.Fs = 44100;
        end
    end
    
end

