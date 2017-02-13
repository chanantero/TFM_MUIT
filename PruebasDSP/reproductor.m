classdef reproductor < matlab.System
     
    properties(Nontunable, SetAccess = private)
        frameSizeReading
        frameSizeWriting
        driver
        device
    end
    
    properties
        getDelayFun
        getAttenFun
    end
    
    properties(DiscreteState)
    end
    
    properties(SetAccess = private, SetObservable, AbortSet)
        playingState
        numChannels % Current number of output channels. It will be set each time the writing device changes
    end

    properties(SetAccess = private)
        audioFileName     
        fileReader % This propertie should be private
    end

    properties(Access = private)

        processor
        player      
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
        
        function updateNumOutputChannels(obj)
            % The number of channels is the maximum number of
            % channels of the writing device
            aux = audioDeviceWriter('Device', obj.device);
            inf = info(aux);
            numChann = inf.MaximumOutputChannels;
            
            obj.numChannels = numChann;
        end
        
        function setProps(obj, frameSize, device, driver)
            if obj.playingState == playingStateClass('stopped')
                obj.frameSizeReading = frameSize;
                obj.frameSizeWriting = frameSize;
                obj.device = device;
                obj.driver = driver;
            end
        end
        
        function stop(obj)
            release(obj.fileReader);
            release(obj.processor);
            release(obj.player);
            obj.playingState = playingStateClass('stopped');
        end
        
        function setDefaultProperties(obj)
            obj.frameSizeReading = 1024*10;
            obj.frameSizeWriting = 1024*10;
            obj.driver = 'DirectSound';
            obj.device = 'Default';
        end
        
        function reproduce(obj)
            % Timing control
            margin = 0.01;
            counter = 0;
            minBufferDepletionTime = 0;
            t0 = tic;
            
            while strcmp(obj.playingState, 'playing') && ~isDone(obj.fileReader) % Only reproduce if it is playing
                t = toc(t0);
                
                if t >= minBufferDepletionTime
                    delay = obj.getDelayFun();
                    step(obj, delay);
                    
                    % Timing control
                    counter = counter + 1;
                    minBufferDepletionTime = counter*obj.frameSizeReading/obj.player.Fs - margin;
                    t = toc(t0);
                    pause(minBufferDepletionTime - t);
                else
                    pause(0.01)
                end
            end
        end
        
    end
        
    methods(Access = public)
        
        function obj = reproductor(fig, position)
            obj.propPanel = propertiesPanel(obj, fig, position, obj.setFrameSize, obj.setDevice, obj.setDriver);
            
            obj.setDefaultProperties();
            obj.playingState = playingStateClass('stopped');
            obj.audioFileName = '';
            
            obj.getDelayFun = @() obj.getDelay();
            
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
            if isempty(order);   return;     end
            
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
                            obj.playingState = playingStateClass('playing');
                            obj.reproduce();
                            % Then, it will play again sincee the state
                            % hasn't been changed
                        case 'resume'
                            % Continue song
                            obj.playingState = playingStateClass('playing');
                            obj.reproduce();
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
        
        function setDriver(obj, driver)
            % Get devices for the new driver
            aux = audioDeviceWriter('Driver', driver);
            audioDevices = getAudioDevices(aux);
            
            obj.setProps(obj.frameSizeReading, audioDevices{1}, driver);
            obj.updateNumOutputChannels();
        end
        
        function setDevice(obj, device)
            obj.setProps(obj.frameSizeReading, device, obj.driver);
            obj.updateNumOutputChannels();
        end
        
        function setFrameSize(obj, frameSize)
            obj.setProps(frameSize, obj.device, obj.driver);
        end
        
        function devices = getWritingDevices(obj)
            aux = audioDeviceWriter('Driver', obj.driver);
            devices = getAudioDevices(aux);
            % devices = obj.player.getAvailableDevices();
        end
        
    end
    
end

