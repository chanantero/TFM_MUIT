classdef reproductorRecorder < matlab.System
    
    % The user can set this properties when the object is not locked, i.e,
    % when playingState is stopped
    properties(Nontunable, SetAccess = private)
        frameDuration % In seconds
        comMatrix % Commutation matrix. numReaders x numPlayers.
        
        % Signal provider settings
        mode % file, sinusoidalm customSignal, func
        FsGenerator % Sampling frequency (For mode == sinusoidal or mode == customSignal or mode == func)
        % mode == file
        audioFileName % Cell string array with as many elements as signalReaders
        % mode == sinusoidal
        amplitude
        phase
        frequency
        % mode == customSignal
        customSignal
        % mode == func
        signalFunc
        
        % Process settings
        enableProc % Enable processing
        modeProc % realTime or predefined
        % modeProc == realTime
        getDelayFun % Cell array of functions with as many elements as processors
        getAttenFun % Cell array of functions with as many elements as processors
        % modeProc == predefined
        pred_t
        pred_delay
        pred_atten
        
        % Player settings
        driver % Cell string array with as many elements as players
        device % Cell string array with as many elements as players
        Fs_player
        defaultChannelMapping
        channelMapping_player
        
        % Recorder settings
        driver_recorder % Cell string array with as many elements as recorders
        device_recorder % Cell string array with as many elements as recorders
        Fs_recorder
        channelMapping_recorder
    end
    
    % Tunable properties
    properties(SetAccess = private)
        comMatrixCoef % Coefficients of the commutation matrix
    end
    
    % Derivated properties from Non-Tunnable properties that are used to
    % set up the properties of the reader, processor and player objects
    properties(SetAccess = private)
        % Reader settings
        frameSizeReading % Array with as many elements as signalReaders
        
        % Player settings
        frameSizeWriting % Array with as many elements as players
        
        % Recorder settings
        frameSizeRecorder % Array with as many elements as recorders
    end
    
    properties(DiscreteState, SetObservable)
        count
    end
    
    % Reading properties
    % The user cannot set this properties directly, but are useful information that
    % can be viewed by other objects
    properties(SetAccess = private, SetObservable, AbortSet)
        playingState % Playing state
        numChannels % Number of channels of the player objects
    end
    
    properties(Dependent, SetAccess = private)
        Fs_reader % Sampling frequency of readers.
        numReaders
        numPlayers
        numLinks
        numRecorders 
        numRecorderChannels % Number of channels of the recorder objects
        recorded
    end
    
    % Private properties: objects
    properties(SetAccess = private)
        signalReader % Cell Array
        procParamProvider % Cell Array
        processor % Cell Array as non-null elements of comMatrix
        player % Cell Array
        recorder % Cell Array
    end
    
    % System object methods
    methods(Access = protected)
        
        function setupImpl(obj)
                    
            obj.setPropertiesSignalProviders();
            
            obj.setPropertiesPlayers();
            
            obj.setPropertiesProcessors();
            
            obj.setPropertiesRecorder();
                   
            for k = 1:obj.numReaders
                setup(obj.signalReader{k})
            end
            
            for k = 1:obj.numLinks
                setup(obj.processor{k}, [], [], [])
                setup(obj.procParamProvider{k});
            end
            
            for k = 1:obj.numPlayers
                setup(obj.player{k}, [])
            end
            
            for k = 1:obj.numRecorders
                setup(obj.recorder{k})
            end
        end
        
        function [numUnderrun, tics] = stepImpl(obj, delays, attenuations)
            
            % Useful parameters
            numReader = obj.numReaders;
            numPlayer = obj.numPlayers;
            indCom = obj.getLinkAbsInd();
            [indReader, indPlayer] = obj.getLinkSubInd();
            numLink = obj.numLinks;
            numRec = obj.numRecorders;
                       
            % Read from file
            audioInput = cell(numReader, 1);
            for k = 1:numReader
                audioInput{k} = step(obj.signalReader{k});
                if obj.enableProc
                    audioInput{k} = mean(audioInput{k}, 2); % From Stereo to Mono
                end
            end
            
            % Process
            audioProc = cell(numLink, 1);
            for k = 1:numLink
                if obj.enableProc
                    audioProc{k} = step(obj.processor{k}, audioInput{indReader(k)}, delays{k}, attenuations{k});
                else
                    audioProc{k} = audioInput{indReader(k)};
                end
            end
            
            % Transform Frequency
            audioProcFs = cell(numLink, 1);
            for k = 1:numLink
                audioProcFs{k} = resample(audioProc{k}, obj.frameSizeWriting(indPlayer(k)), obj.frameSizeReading(indReader(k)));
            end
            
            % Sum outputs according to the commutation matrix
            % All processed signals corresponding to the same player must
            % have the same number of samples
            audioOutput = cell(numPlayer, 1);
            for k = 1:numPlayer
                ind = find(indPlayer == k);
                audioOutput{k} = zeros(obj.frameSizeWriting(k), obj.numChannels(k));
                for l = 1:numel(ind)
                    audioOutput{k} = audioOutput{k} + obj.comMatrixCoef(indCom(ind(l)))*audioProcFs{ind(l)};
                end
            end
            
            % Write to audio device buffer
            indPlay = obj.getActivePlayers();
            numUnderrun = zeros(numPlayer, 1);
            tics = uint64(zeros(numPlayer, 1)); % When does it finish processing?
            for k = indPlay
                [numUnderrun(k), tics(k)] = step(obj.player{k}, audioOutput{k});
            end
                         
            % Record
            for k = 1:numRec
                step(obj.recorder{k});
            end
            
            % Update discrete state
            obj.count = obj.count + 1;
        end
        
        function resetImpl(obj)
            obj.reset_priv();
        end
        
        function releaseImpl(obj)
            obj.stop();
        end
        
        function validatePropertiesImpl(obj)
            for k = 1:obj.numReaders
                if obj.mode(k) == originType('file')
                    % Check that the audioFileName exists
                    a = dir(obj.audioFileName{k});
                    assert(numel(a) > 0, 'reproductor:wrongProperty', 'The property audioFileName must specify an audio file that exists')
                    
                    % Check that the audioFileName has the right extension
                    [~, ~, ext] = fileparts(obj.audioFileName{k});
                    extensions = {'.wav', '.mp3'};
                    assert(ismember(ext, extensions), 'reproductor:wrongProperty', ['The property audioFileName must have one of the next extensions: ', strjoin(extensions, ', ')])
                end
            end
        end
        
    end
    
    % Getters and setters
    methods
        
        function numPlayers = get.numPlayers(obj)
            numPlayers = numel(obj.player);
        end
        
        function numReaders = get.numReaders(obj)
            numReaders = numel(obj.signalReader);
        end
        
        function numLinks = get.numLinks(obj)
            numLinks = sum(obj.comMatrix(:));
        end
        
        function numRecorders = get.numRecorders(obj)
            numRecorders = numel(obj.recorder);
        end
        
        function numRecorderChannels = get.numRecorderChannels(obj)
            N = obj.numRecorders;
            numRecorderChannels = zeros(N, 1);
            for k = 1:N
                numRecorderChannels(k) = obj.recorder{k}.NumChannels;
            end
           
        end
        
        function recorded = get.recorded(obj)
            N = obj.numRecorders;
            recorded = cell(N, 1);
            for k = 1:N
                recorded{k} = obj.recorder{k}.recorded;
            end
        end
        
        function Fs = get.Fs_reader(obj)
            Fs = zeros(obj.numReaders, 1);
            for k = 1:obj.numReaders
                Fs(k) = obj.signalReader{k}.SampleRate;
            end
        end
        
    end
    
    
    methods(Access = private)
        
        function setPropertiesSignalProviders(obj, index)
            if nargin == 1
                index = 1:obj.numReaders;
            end
            
            for k = index
                obj.signalReader{k}.FileName = obj.audioFileName{k};
                obj.signalReader{k}.mode = obj.mode(k);
                obj.signalReader{k}.amplitude = obj.amplitude(k);
                obj.signalReader{k}.phase = obj.phase(k);
                obj.signalReader{k}.frequency = obj.frequency(k);
                obj.signalReader{k}.SampleRate = obj.FsGenerator(k);
                obj.signalReader{k}.customSignal = obj.customSignal{k};
                obj.signalReader{k}.signalFunc = obj.signalFunc{k};
            end
            
            obj.frameSizeReading = floor(obj.frameDuration*obj.Fs_reader);
            
            for k = index
                obj.signalReader{k}.SamplesPerFrame = obj.frameSizeReading(k);
            end
        end
        
        function setPropertiesProcessors(obj, index)
            if nargin == 1
                index = 1:obj.numLinks;
            end
            
            [indReader, indPlayer] = obj.linkAbs2Sub(index);
            
            % Set properties for the processors and processor parameter
            % providers
            for k = 1:numel(index)
                obj.processor{k}.Fs = obj.Fs_reader(indReader(k));
                obj.processor{k}.numChannels = obj.numChannels(indPlayer(k));
                
                obj.procParamProvider{k}.mode = obj.modeProc(k);
                obj.procParamProvider{k}.FrameSize = obj.frameSizeReading(indReader(k));
                
                obj.procParamProvider{k}.getDelayFun = obj.getDelayFun{k};
                obj.procParamProvider{k}.getAttenFun = obj.getAttenFun{k};
                
                obj.procParamProvider{k}.SampleRate = obj.Fs_reader(indReader(k));
                obj.procParamProvider{k}.t_change = obj.pred_t{k};
                obj.procParamProvider{k}.delays = obj.pred_delay{k};
                obj.procParamProvider{k}.attenuations = obj.pred_atten{k};
            end
      
        end
        
        function setPropertiesPlayers(obj, index)
            if nargin == 1
                index = 1:obj.numPlayers;
            end

            obj.frameSizeWriting = floor(obj.frameDuration*obj.Fs_player);
            
            for k = index
                obj.player{k}.driver = obj.driver{k};
                obj.player{k}.device = obj.device{k};
                obj.player{k}.frameSize = obj.frameSizeWriting(k);
                obj.player{k}.Fs = obj.Fs_player(k);
                obj.player{k}.DefaultChannels = obj.defaultChannelMapping(k);
                obj.player{k}.channelMapping = obj.channelMapping_player{k};
            end
        end
      
        function setPropertiesRecorder(obj, index)
             if nargin == 1
                index = 1:obj.numRecorders;
             end
             
             obj.frameSizeRecorder = floor(obj.frameDuration*obj.Fs_recorder);
            
             for k = index
                 obj.recorder{k}.Driver = obj.driver_recorder{k};
                 obj.recorder{k}.Device = obj.device_recorder{k};
                 obj.recorder{k}.SamplesPerFrame = obj.frameSizeRecorder(k);
                 obj.recorder{k}.SampleRate = obj.Fs_recorder(k);
                 obj.recorder{k}.channelMapping = obj.channelMapping_recorder{k};
             end
        end
        
        function updateNumOutputChannels(obj, index)
           
            if nargin == 1
                index = 1:obj.numPlayers;
            end
            
            numChann = zeros(numel(index), 1);
            for k = 1:numel(index)
                numChann(k) = obj.player{k}.numChannels;
            end
            
            obj.numChannels(index) = numChann;
        end
        
        function stop(obj)
            for k = 1:obj.numReaders
                release(obj.signalReader{k})
            end
            
            for k = 1:obj.numLinks
                release(obj.processor{k})
                release(obj.procParamProvider{k})
            end
            
            for k = 1:obj.numPlayers
                release(obj.player{k})
            end
            
            for k = 1:obj.numRecorders
                release(obj.recorder{k})
            end
            
            obj.playingState = playingStateClass('stopped');
        end
        
        function reset_priv(obj)
            for k = 1:obj.numReaders
                reset(obj.signalReader{k})
            end
            
            for k = 1:obj.numLinks
                reset(obj.processor{k})
                reset(obj.procParamProvider{k});
            end
            
            for k = 1:obj.numPlayers
                reset(obj.player{k})
            end
            
            for k = 1:obj.numRecorders
                reset(obj.recorder{k})
            end
            
            obj.count = 0;
        end
        
        function setReproductorDefaultProperties(obj)
            
            for k = 1:obj.numReaders
                obj.setProps('audioFileName', '', k);
                obj.setProps('frameSizeReading', 1024*10, k);
                obj.setProps('mode', originType('file'), k);
                obj.setProps('amplitude', 1, k);
                obj.setProps('phase', 0, k);
                obj.setProps('frequency', 1, k);
                obj.setProps('FsGenerator', 44100, k);
                obj.setProps('customSignal', [], k);
            end
            
            for k = 1:obj.numPlayers
                obj.setProps('Fs_player', 44100, k);
                obj.setProps('driver', 'DirectSound', k);
                obj.setProps('device', 'Default', k);
                obj.setProps('defaultChannelMapping', true);
            end
                        
            [readInd, playInd] = obj.getLinkSubInd();
            for k = 1:obj.numLinks
                ind = [readInd(k), playInd(k)];
                obj.setProps('modeProc', timeInteractionTypes('realTime'), ind);
                obj.setProps('getDelayFun', [], ind);
                obj.setProps('getAttenFun', [], ind);
                obj.setProps('pred_t', [], ind);
                obj.setProps('pred_delay', [], ind);
                obj.setProps('pred_atten', [], ind);
            end
            
        end
 
        function setRecorderDefaultProperties(obj)         
            for k = 1:obj.numRecorders
                obj.setProps('driver_recorder', 'DirectSound', k);
                obj.setProps('device_recorder', 'Default', k);
                obj.setProps('Fs_recorder', 44100, k);
                obj.setProps('channelMapping_recorder', 1, k);
            end
        end
        
        function setOtherDefaultProperties(obj)
            obj.setProps('frameDuration', 0.25);
        end
        
        function setCommutationMatrixCoef(obj, comMatCoef)
            obj.comMatrixCoef = comMatCoef;
        end
        
        function setCommutationMatrix(obj, comMat)
            obj.comMatrix = comMat;
            
            [numReaders_new, numPlayers_new] = size(comMat);
            numLinks_new = sum(comMat(:));

            obj.comMatrixCoef = double(comMat);
                        
            obj.signalReader = cell(numReaders_new, 1);
            for k = 1:numReaders_new
                obj.signalReader{k} = signalProvider();
            end
            
            obj.processor = cell(numLinks_new, 1);
            obj.procParamProvider = cell(numLinks_new, 1);
            for k = 1:numLinks_new
                obj.procParamProvider{k} = delayAndAttenProvider();
                obj.processor{k} = processSignal();
            end
            
            obj.player = cell(numPlayers_new, 1);
            for k = 1:numPlayers_new
                obj.player{k} = audioPlayer;
            end
            
            obj.setReproductorDefaultProperties();
        end
        
        function setNumRecorders(obj, numRecorders)
            obj.recorder = cell(numRecorders, 1);
            for k = 1:numRecorders
                obj.recorder{k} = audioRecorder;
            end
            
            obj.setRecorderDefaultProperties();
        end
        
        function t_br = delayBetweenDevices(obj, ts_bufferQueueLoad, numUnderruns)
            % Assume that the buffer size and the size of each load to the
            % buffer queue are equal, i.e. the VariableInputSize of the
            % deviceWriterObject is set to false.
            N = numel(ts_bufferQueueLoad); % Number of devices
            
            t_br = zeros(N, 1); % Beginning time of reproduction
            for k = 1:N
                numUnderrunFrames = numUnderruns{k}/obj.frameSizeWriting(k);
                delayLimits = processBufferResults( ts_bufferQueueLoad{k}, numUnderrunFrames, obj.frameDuration, obj.frameDuration );
                delay = delayLimits(1); % Let's take an approximation of the real delay
                t0 = ts_bufferQueueLoad{k}(1); % Reference time against which the delay Limits are calculated
                t_br(k) = delay + t0;
            end
            
            % Normalize
            t_br = t_br - min(t_br);
            
        end
        
        function reproduce(obj)
            if ~strcmp(obj.playingState, 'playing'),  return;  end
            
            numPlay = obj.numPlayers;
            
            % Timing control
            margin = 0.5;
            counter = 0;
            minPause = 0.01;
            minBufferDepletionTime = 0;
            t0 = tic; tref = t0;
            % Buffer control
            offset = zeros(numPlay, 1);
            t_eps = cell(numPlay, 1); % Times of the ending of loading to the buffer queue
            numUnderruns = cell(numPlay, 1);
            % Recorder timing
            del = 2; % Delay of the reproduction with respect to the recording in seconds
                        
            [~, playerIndex] = obj.getLinkSubInd();
            
            for k = 1:obj.numRecorders
                step(obj.recorder{k});
            end
            
            finish = false;
            while ~finish % Only reproduce if it is playing
                t = tic;
                
                if toc(t0) >= minBufferDepletionTime
                    
                    % Get processor parameters: delays and attenuations
                    delays = cell(obj.numLinks, 1);
                    attenuations = cell(obj.numLinks, 1);
                    if obj.enableProc
                        for k = 1:obj.numLinks
                            [delay, attenuations{k}] = step(obj.procParamProvider{k});
                            delay = delay + offset(playerIndex(k));
                            delays{k} = delay;
                        end
                    end
                    
                    
                    % Step
                    try
                        [numUnderrun, t_ep] = step(obj, delays, attenuations);
                                                
                        % Delay between writing devices control
                        for k = 1:obj.numPlayers
                            if numel(t_eps{k}) == 20
                                t_eps{k}(1:end-1) = t_eps{k}(2:end); % Shift
                                t_eps{k}(end) = toc(tref) - toc(t_ep(k)); % New value to the end

                                numUnderruns{k}(1:end-1) = numUnderruns{k}(2:end); % Shift
                                numUnderruns{k}(end) = numUnderrun(k); % New value to the end
                            else
                                t_eps{k} = [t_eps{k}; toc(tref) - toc(t_ep(k))];
                                numUnderruns{k} = [numUnderruns{k}; numUnderrun(k)];
                            end
                        end
                        
                        if numel(t_eps{1}) > 1
                            offset = obj.delayBetweenDevices(t_eps, numUnderruns);
                            offset = max(offset) - offset;
                        end
                    catch e
                        disp(e)
                        warning('There was some error with the step function of reproductor')
                        release(obj);
                        return;
                    end
                    
                    % Timing control
                    if any(numUnderrun > 0)
                        % Interruption in the reproduction. Reset timer
                        fprintf('Underrun. count = %d', obj.count);
                        t0 = t;
                        counter = 0;
                    end
                    
                    counter = counter + 1;
                    minBufferDepletionTime = counter*obj.frameDuration - margin;
                    pause(max(minPause, minBufferDepletionTime - toc(t0)));
                    %                     fprintf(max(minPause, minBufferDepletionTime - toc(t0)))
                else
                    pause(minPause);
                end
                
                % Finish condition
                if ~strcmp(obj.playingState, 'playing')
                    finish = true;
                else
                    done = false(obj.numReaders, 1);
                    for k = 1:obj.numReaders
                        done(k) = isDone(obj.signalReader{k});
                    end
                    if all(done)
               
                        % Record to overcompensate delays
                        N = ceil(del/obj.frameDuration);
                        for k = 1:obj.numRecorders
                            for n = 1:N
                                step(obj.recorder{k});
                            end
                        end
                        
                        finish = true;
                        release(obj)
                    end
                end
                
            end
            
        end
        
        function lowLevelOrder(obj, action)
                      
            % There are 5 basic actions:
            % - lock
            % - unlock
            % - resume
            % - pause
            % - restart
            
            switch obj.playingState
                case playingStateClass('playing') % Reproducing (obviously locked)
                    switch action
                        case 'lock'
                            % It is already locked, so do nothing
                        case 'unlock'
                            assert(isLocked(obj), 'reproductor:WrongState', 'When the state is playing, the object should be locked');
                            release(obj); % It only executes if obj is locked. Make sure that always the state is playing, the obect is locked
                        case 'resume'
                            % It is already reproducing, so do nothing
                        case 'pause'
                            % Change state. The rest is left the same way
                            obj.playingState = playingStateClass('paused');
                        case 'restart'
                            reset(obj);
                        case ''
                            % Null action
                    end
                    
                case playingStateClass('stopped') % Unlocked
                    switch action
                        case 'lock'
                            setup(obj, {}, {});
                            reset(obj);
                            obj.playingState = playingStateClass('paused');
                        case 'unlock'
                            % It is already stopped, so do nothing
                        case 'resume'
                            % The only possible action is locking, so do
                            % nothing
                        case 'pause'
                            % The only possible action is locking, so do
                            % nothing
                        case 'restart'
                            % The only possible action is locking, so do
                            % nothing
                        case ''
                            % Null action
                    end
                    
                case playingStateClass('paused') % Locked but not reproducing
                    switch action
                        case 'lock'
                            % It is already locked, so do nothing
                        case 'unlock'
                            assert(isLocked(obj), 'reproductor:WrongState', 'When the state is playing, the object should be locked');
                            release(obj); % It only executes if obj is locked. Make sure that always the state is playing, the obect is locked
                        case 'resume'
                            % Continue song
                            obj.playingState = playingStateClass('playing');
                            obj.reproduce();
                        case 'pause'
                            % It is already paused, so do nothing
                        case 'restart'
                            reset(obj);
                        case ''
                            % Null action
                    end
                    
                otherwise
                    error('reproductor:nonValidPlayingState', 'The playing state is not valid')
            end
        end
        
    end
    
    
    methods(Access = public)
        
        function obj = reproductorRecorder()
            
            obj.playingState = playingStateClass('stopped');
            
            % Reading object
            obj.signalReader = {signalProvider()};
            
            % Processing object
            obj.enableProc = true;
            obj.processor = {processSignal()};
            obj.procParamProvider = {delayAndAttenProvider()};
            
            % Writing object
            obj.player = {audioPlayer};
            
            obj.comMatrix = true;
            obj.comMatrixCoef = 1;
            obj.setReproductorDefaultProperties();
            
            % Recording object
            obj.setNumRecorders(1);
            
            obj.setOtherDefaultProperties();
            
        end
        
        function executeOrder(obj, action)
                        
            switch obj.playingState
                case playingStateClass('playing')
                    switch action
                        case 'stop'
                            obj.lowLevelOrder('unlock');
                        case 'play'
                            % Play from the beggining
                            obj.lowLevelOrder('restart');   
                        case 'resume'
                            % Do nothing
                        case 'pause'
                            % Change state. The rest is left the same way
                            obj.lowLevelOrder('pause');
                        case ''
                            % Null action
                    end
                    
                case playingStateClass('stopped')
                    switch action
                        case 'play'
                            % Play from the beginning
                            obj.lowLevelOrder('lock');
                            obj.lowLevelOrder('resume');
                        case 'resume'
                            % Do nothing
                        case 'stop'
                            % Do nothing
                        case 'pause'
                            % Do nothing
                        case ''
                            % Null action
                    end
                    
                case playingStateClass('paused')
                    switch action
                        case 'play'
                            % Start current song from the begginging
                            obj.lowLevelOrder('restart');
                            obj.lowLevelOrder('resume');
                        case 'resume'
                            % Continue song
                            obj.lowLevelOrder('resume');
                        case 'pause'
                            % Do nothing
                        case 'stop'
                            % Completely release
                            obj.lowLevelOrder('unlock');
                        case ''
                            % Null action
                    end
                    
                otherwise
                    error('reproductor:nonValidPlayingState', 'The playing state is not valid')
            end
            
        end
        
        function setProps(obj, parameters, values, indices )
            
            % Parse inputs
            if ~iscell(parameters)
                parameters = {parameters};
                values = {values};
                
                if nargin < 4
                    indices = {0};
                else
                    indices = {indices};
                end
            end
            
            if obj.playingState == playingStateClass('stopped')
        
                % Set parameters
                for k = 1:numel(parameters)
                    
                    parameter = parameters{k};
                    value = values{k};
                    index = indices{k};
                    
                    switch parameter
                        case 'audioFileName'
                            obj.(parameter){index} = value;
                        case 'driver'
                            obj.(parameter){index} = value;
                            
                            % Set the value on the correspondent player
                            player_this = obj.player{index};
                            player_this.driver = value;
                                                        
                            % Get devices for the new driver and set it
                            audioDevices = player_this.getAudioDevices();
                            obj.setProps('device', audioDevices{1}, index);
                            
                        case 'device'
                            obj.(parameter){index} = value;
                            
                            % Set the value on the correspondent player and do
                            % the setup and release
                            player_this = obj.player{index};
                            player_this.device = value;
                            setup(player_this, []);
                            release(player_this);
                            
                            % Update the number of channels for that device
                            obj.updateNumOutputChannels(index);
                        case 'modeProc'
                            ind = obj.sub2LinkAbs(index(1), index(2));
                            obj.(parameter)(ind) = value;
                        case 'getDelayFun'                            
                            ind = obj.sub2LinkAbs(index(1), index(2));
                            obj.(parameter){ind} = value;
                        case 'getAttenFun'
                            ind = obj.sub2LinkAbs(index(1), index(2));
                            obj.(parameter){ind} = value;
                        case 'pred_t'
                            ind = obj.sub2LinkAbs(index(1), index(2));
                            obj.(parameter){ind} = value;
                        case 'pred_delay'
                            ind = obj.sub2LinkAbs(index(1), index(2));
                            obj.(parameter){ind} = value;
                        case 'pred_atten'
                            ind = obj.sub2LinkAbs(index(1), index(2));
                            obj.(parameter){ind} = value;
                        case 'comMatrix'
                            obj.setCommutationMatrix(value);
                        case 'comMatrixCoef'
                            obj.setCommutationMatrixCoef(value);
                        case 'frameSizeReading'
                            obj.(parameter)(index) = value;
                        case 'Fs_player'
                            obj.(parameter)(index) = value;
                        case 'frameDuration'
                            obj.(parameter) = value;
                        case 'mode'
                            obj.(parameter)(index) = value;
                        case 'amplitude'
                            obj.(parameter)(index) = value;
                        case 'phase'
                            obj.(parameter)(index) = value;
                        case 'frequency'
                            obj.(parameter)(index) = value;
                        case 'FsGenerator'
                            obj.(parameter)(index) = value;
                        case 'numRecorders'
                            obj.setNumRecorders(value);
                        case 'driver_recorder'
                            obj.(parameter){index} = value;
                            
                            % Set the value on the correspondent recorder
                            recorder_this = obj.recorder{index};
                            recorder_this.Driver = value;
                                                        
                            % Get devices for the new driver and set it
                            audioDevices = recorder_this.getAudioDevices();
                            obj.setProps('device_recorder', audioDevices{1}, index);
                            
                        case 'device_recorder'
                            obj.(parameter){index} = value;
                        case 'Fs_recorder'
                            obj.(parameter)(index) = value;
                        case 'customSignal'
                            obj.(parameter){index} = value;
                        case 'signalFunc'
                            obj.(parameter){index} = value;
                        case 'enableProc'
                            obj.(parameter) = value;
                        case 'channelMapping_recorder'
                            obj.(parameter){index} = value;
                        case 'defaultChannelMapping'
                            obj.(parameter)(index) = value;
                        case 'channelMapping_player'
                            obj.(parameter){index} = value;
                        otherwise
                            error('reproductor_plus:setProps', 'The first argument must be the name of an existing parameter')
                    end
                    
                end
            else
                % You can modify tunable properties
                for k = 1:numel(parameters)
                    parameter = parameters{k};
                    value = values{k};
                    
                    switch parameter
                        case 'comMatrixCoef'
                            obj.setCommutationMatrixCoef(value);
                        otherwise
                            error('reproductor_plus:setProps', 'The first argument must be the name of an existing parameter')
                    end
                end
            end
        end
          
        function indices = getLinkAbsInd(obj)
            indices = find(obj.comMatrix);
        end
        
        function [readerIndex, playerIndex] = getLinkSubInd(obj)
            absLinkInd = obj.getLinkAbsInd();
            [readerIndex, playerIndex] = ind2sub(size(obj.comMatrix), absLinkInd);
        end
        
        function indices = sub2LinkAbs(obj, readerIndex, playerIndex)
            absInd = sub2ind(size(obj.comMatrix), readerIndex, playerIndex);
            absLinkInd = obj.getLinkAbsInd();
            [~, indices] = ismember(absInd, absLinkInd);
        end
             
        function [readerIndex, playerIndex] = linkAbs2Sub(obj, index)
            absInd = obj.getLinkAbsInd();
            [readerIndex, playerIndex] = ind2sub(size(obj.comMatrix), absInd(index));
        end
        
        function indices = getActivePlayers(obj)
            indices = find(any(obj.comMatrix, 1));
        end
        
        function indices = getActiveReaders(obj)
            indices = find(any(obj.comMatrix, 2));
        end
    end
    
end

