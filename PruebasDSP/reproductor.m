classdef reproductor < matlab.System
     
    properties(Nontunable)
        frameSizeReading
        frameSizeWriting
        driver
        device
    end
    
    properties
        givenControl
    end
    
    properties(DiscreteState)
    end

    properties(SetAccess = private)
        playingState
        audioFileName            
    end

    properties(Access = private)
        fileReader
        processor
        player
        numChannels % Current number of output channels. It will be set during the setup
    end
    
    methods(Access = protected)
        
        function setupImpl(obj)
            obj.fileReader.Filename = obj.audioFileName;
            obj.fileReader.SamplesPerFrame = obj.frameSizeReading;
            
            obj.player.frameSize = obj.frameSizeWriting;
            obj.player.driver = obj.driver;
            obj.player.device = obj.device;
            
            % The sampling frequency will depend on the sampling frequency
            % of the original file
            Fs = obj.fileReader.SampleRate;
            obj.player.Fs = Fs;
            obj.processor.Fs = Fs;
            
            % The number of channels will depend on the maximum number of
            % channels of the writing device
            aux = audioDeviceWriter('Device', obj.device);
            inf = info(aux);
            numChann = inf.MaximumOutputChannels;
            obj.numChannels = numChann;
            obj.processor.numChannels = numChann;
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
            reset(obj.fileReader);
            reset(obj.processor);
            reset(obj.player);
        end

        function releaseImpl(obj)
            obj.stop();
        end
        
        
    end
    
    methods(Access = private)
        function delays = getDelay(obj)
            delays = zeros(obj.numChannels, 1);
        end
    end
    
    methods(Access = public)
        
        function reproduce(obj)
            while ~isDone(obj.fileReader)
                % Only reproduce if it is playing
                if strcmp(obj.playingState, 'playing')
                    delay = obj.getDelay();
                    step(obj, delay);
                    pause(0.01)
                else
                    return;
                end
            end
        end
        
        function obj = reproductor()
            obj.setDefaultProperties();
            obj.playingState = playingStateClass('stopped');
            obj.givenControl = false;
            
            % Reading object
            obj.fileReader = dsp.AudioFileReader();
            
            % Processing object
            obj.processor = processSignal('variable', true, 'delayType', 'forward');
            
            % Writing object
            obj.player = audioPlayer;
            
        end
                
        function executeOrder(obj, order)
            % 'order' is a structure with 1 or 2 fields.
            % The field 'action' is mandatory
            % The field 'fileName' is not mandatory. It refers to the
            % active track
            p = inputParser;
            addParameter(p, 'action', '');
            addParameter(p, 'fileName', '');
            parse(p, order);
            
            action = p.Results.action;
            
            switch obj.playingState
                case playingStateClass('playing')
                    switch action
                        case 'stop'
                            release(obj);
                        case 'play'
                            % Two possibilities. Play the same song (reset)
                            % or play a new song
                            if ismember('fileName', p.UsingDefaults)
                                % Play current song from the beggining
                                % First, stop but without allowing the tuning
                                % of anything, i.e., don't release the objects
                                reset(obj);
                                % Then, it will play again sincee the state
                                % hasn't been changed
                            else
                                % Play a new audio file
                                % First, stop
                                release(obj);
                                % Then, assign the new file name
                                obj.audioFileName = p.Results.fileName;
                                % Then, play again
                                obj.playingState = playingStateClass('playing');
                                setup(obj, []);
                            end
                        case 'resume'
                            % Do nothing
                        case 'pause'
                            % Change state. The rest is left the same way
                            obj.playingState = playingStateClass('pause');
                        case 'assignTrack'
                            % Play a new audio file
                            % First, stop
                            release(obj);
                            % Then, assign the new file name
                            obj.audioFileName = p.Results.fileName;
                            % Then, play again
                            obj.playingState = playingStateClass('playing');
                            setup(obj, []);
                        case ''
                            % Null action
                    end
                    
                case playingStateClass('stopped')
                    switch action
                        case 'play'
                            % Two possibilities. Play the same song (reset)
                            % or play a new song
                            if ~ismember('fileName', p.UsingDefaults)
                                % Play a new audio file
                                % Assign the new file name
                                obj.audioFileName = p.Results.fileName;
                            end
                            setup(obj, []);
                            obj.playingState = playingStateClass('playing');
                            obj.reproduce();
                        case 'resume'
                            % The same as play if there was no song
                            % specified
                            setup(obj, []);
                            obj.playingState = playingStateClass('playing');
                            obj.reproduce();
                        case 'stop'
                            % Do nothing
                        case 'pause'
                            % Do nothing
                        case 'assignTrack'
                            % Assign the new file name
                            obj.audioFileName = p.Results.fileName;
                        case ''
                            % Null action
                    end
                    
                case playingStateClass('paused')
                    switch action
                        case 'play'
                            % Start current song from the begginging
                            % First, stop but without allowing the tuning
                            % of anything, i.e., don't release the objects,
                            % only reset it
                            reset(obj);
                            obj.playingState = playingStateclass('playing');
                            obj.reproduce();
                            % Then, it will play again sincee the state
                            % hasn't been changed
                        case 'resume'
                            % Continue song
                            obj.playingState = playingStateClass('playing');
                        case 'pause'
                            % Do nothing
                        case 'stop'
                            % Completely release
                            release(obj);
                        case 'assignTrack'
                            % New audio file, but don't start reproducing
                            % First, stop
                            release(obj);
                            % Then, assign the new file name
                            obj.audioFileName = p.Results.fileName;
                            obj.playingState = playingStateClass('stopped');
                        case ''
                            % Null action
                    end
                
                otherwise
                    error('reproductor:nonValidPlayingState', 'The playing state is not valid')
            end
                    
        end      
                
        function stop(obj)
            release(obj.fileReader);
            release(obj.processor);
            release(obj.player);
            obj.playingState = playingStateClass('stopped');
        end
        
        function setDefaultProperties(obj)
            obj.frameSizeReading = 1024;
            obj.frameSizeWriting = 1024;
            obj.driver = 'DirectSound';
            obj.device = 'Default';
        end
    end
    
end

