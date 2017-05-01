classdef reproductor_plus < matlab.System
    
    % The user can set this properties when the object is not locked, i.e,
    % when playingState is stopped
    properties(Nontunable, SetAccess = private)
        frameDuration % In seconds
        
        % Reader settings
        audioFileName % Cell string array with as many elements as fileReaders
        
        % Processor settings
        getDelayFun % Cell array of functions with as many elements as players
        getAttenFun % Cell array of functions with as many elements as players
        comMatrix % Commutation matrix. numReaders x numPlayers.
        
        % Player settings
        driver % Cell string array with as many elements as players
        device % Cell string array with as many elements as players
        Fs_player
    end
    
    % Derivated properties from Non-Tunnable properties that are used to
    % set up the properties of the reader, processor and player objects
    properties(SetAccess = private)
        % Reader settings
        frameSizeReading % Array with as many elements as fileReaders
        
        % Player settings
        frameSizeWriting % Array with as many elements as players
    end
    
    properties(DiscreteState)
        count
    end
    
    % Reading properties
    % The user cannot set this properties directly, but are useful information that
    % can be viewed by other objects
    properties(SetAccess = private, SetObservable)%, AbortSet)
        playingState % Playing state
        numChannels % Current number of output channels. It will be set each time a writing device changes
        Fs_reader % Sampling frequency of readers.
        numReaders
        numPlayers
        numLinks
    end
    
    % Private properties: objects
    properties(SetAccess = private)
        fileReader % Cell Array
        processor % Cell Array as non-null elements of comMatrix
        player % Cell Array
    end
    
    % System object methods
    methods(Access = protected)
        
        function setupImpl(obj)
            
%             % Calculate Fs_player automatically
%             firstPlayer = zeros(obj.numPlayers, 1);
%             for k = 1:obj.numPlayers
%                 firstPlayer(k) = find(obj.comMatrix(:, k), 1, 'first');
%             end
%             obj.Fs_player(k) = obj.Fs_reader(firstPlayer(k));
            
            obj.setPropertiesFileReaders();
            
            obj.setPropertiesPlayers();
            
            obj.setPropertiesProcessors();
            
            
        end
        
        function [numUnderrun, tics] = stepImpl(obj, delays, attenuations)
            
            % Useful parameters
            numReader = obj.numReaders;
            numPlayer = obj.numPlayers;
            indCom = find(obj.comMatrix ~= 0);
            [indReader, indPlayer] = ind2sub([numReader, numPlayer], indCom);
            numLink = obj.numLinks;
            
            % Read from file
            audioInput = cell(numReader, 1);
            for k = 1:numReader
                audioInput{k} = step(obj.fileReader{k});
                audioInput{k} = mean(audioInput{k}, 2); % From Stereo to Mono
            end
            
            % Process
            audioProc = cell(numLink, 1);
            for k = 1:numLink
                audioProc{k} = step(obj.processor{k}, audioInput{indReader(k)}, delays{k}, attenuations{k});
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
                audioOutput{k} = zeros(size(audioProcFs{ind(1)}));
                for l = 1:numel(ind)
                    audioOutput{k} = audioOutput{k} + obj.comMatrix(indCom(ind(l)))*audioProcFs{ind(l)};
                end
            end
            
            % Write to audio device buffer
            numUnderrun = zeros(numPlayer, 1);
            tics = uint64(zeros(numPlayer, 1)); % When does it finish processing?
            for k = 1:numPlayer
                [numUnderrun(k), tics(k)] = step(obj.player{k}, audioOutput{k});
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
    
    % Getters and setters
    methods
        function numPlayers = get.numPlayers(obj)
            numPlayers = numel(obj.player);
        end
        
        function numReaders = get.numReaders(obj)
            numReaders = numel(obj.fileReader);
        end
        
        function numLinks = get.numLinks(obj)
            numLinks = sum(obj.comMatrix ~= 0);
        end
        
        function Fs = get.Fs_reader(obj)
            Fs = zeros(obj.numReaders, 1);
            for k = 1:obj.numReaders
                Fs(k) = obj.fileReader{k}.SampleRate;
            end
        end
        
    end
    
    
    methods(Access = private)
        
        function setPropertiesFileReaders(obj, index)
            if nargin == 1
                index = 1:obj.numReaders;
            end
            
            for k = index
                obj.fileReader{k}.Filename = obj.audioFileName{k};
            end
            
            obj.frameSizeReading = floor(obj.frameDuration*obj.Fs_reader);
            
            for k = index
                obj.fileReader{k}.SamplesPerFrame = obj.frameSizeReading(k);
            end
        end
        
        function setPropertiesProcessors(obj, index)
            if nargin == 1
                index = 1:obj.numLinks;
            end
            
            [indReader, indPlayer] = obj.linkAbs2Sub(index);
            
            % Set properties for the processors
            for k = 1:numel(index)
                obj.processor{k}.Fs = obj.Fs_reader(indReader(k));
                obj.processor{k}.numChannels = obj.numChannels(indPlayer(k));
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
                release(obj.fileReader{k})
            end
            
            for k = 1:obj.numLinks
                release(obj.processor{k})
            end
            
            for k = 1:obj.numPlayers
                release(obj.player{k})
            end
            
            obj.playingState = playingStateClass('stopped');
        end
        
        function reset_priv(obj)
            for k = 1:obj.numReaders
                reset(obj.fileReader{k})
            end
            
            for k = 1:obj.numLinks
                reset(obj.processor{k})
            end
            
            for k = 1:obj.numPlayers
                reset(obj.player{k})
            end
            
            obj.count = 0;
        end
        
        function setDefaultProperties(obj)
            
            for k = 1:obj.numReaders
                obj.setProps('audioFileName', '', k);
                obj.setProps('frameSizeReading', 1024*10, k);
            end
            
            for k = 1:obj.numPlayers
                obj.setProps('Fs_player', 44100, k);
                obj.setProps('driver', 'DirectSound', k);
                obj.setProps('device', 'Default', k);
            end
            
            obj.setProps('frameDuration', 0.25, k);
            
            [readInd, playInd] = obj.getLinkSubInd();
            for k = 1:obj.numLinks
                ind = [readInd(k), playInd(k)];
                obj.setProps('getDelayFun', [], ind);
                obj.setProps('getAttenFun', [], ind);
            end
            
        end
        
        function setCommutationMatrix(obj, comMat)
            obj.comMatrix = comMat;
            
            [numReaders_new, numPlayers_new] = size(comMat);
            numLinks_new = sum(comMat ~= 0);
            
            obj.fileReader = cell(numReaders_new, 1);
            for k = 1:numReaders_new
                obj.fileReader{k} = dsp.AudioFileReader();
            end
            
            obj.processor = cell(numLinks_new, 1);
            for k = 1:numLinks_new
                obj.processor{k} = processSignal('delayType', 'forward');
            end
            
            obj.player = cell(numPlayers_new, 1);
            for k = 1:numPlayers_new
                obj.player{k} = audioPlayer;
            end
            
            obj.setDefaultProperties();
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
        
        function reproduce(obj)
            
            numPlay = obj.numPlayers;
            
            % Timing control
            margin = 0.01;
            counter = 0;
            minPause = 0.01;
            minBufferDepletionTime = 0;
            t0 = tic; tref = t0;
            % Buffer control
            offset = zeros(numPlay, 1);
            t_eps = cell(numPlay, 1); % Times of the ending of loading to the buffer queue
            numUnderruns = cell(numPlay, 1);
                        
            finish = false;
            while ~finish % Only reproduce if it is playing
                t = tic;
                
                if toc(t0) >= minBufferDepletionTime
                    
                    [readerIndex, playerIndex] = obj.getLinkSubInd();

                    delays = cell(obj.numLinks, 1);
                    attenuations = cell(obj.numLinks, 1);
                    for k = 1:obj.numLinks
                        delay = obj.getDelayFun{k}() + offset(playerIndex(k));
                        attenuation = obj.getAttenFun{k}();
                                          
                        delays{k} = repmat(delay', obj.frameSizeReading(readerIndex(k)), 1);
                        attenuations{k} = repmat(attenuation', obj.frameSizeReading(readerIndex(k)), 1);
                    end
                    
                    try
                        [numUnderrun, t_ep] = step(obj, delays, attenuations);
                        
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
                    catch
                        warning('There was some error with the step function of reproductor')
                        release(obj);
                        return;
                    end
                    
                    % Timing control
                    if any(numUnderrun > 0)
                        % Interruption in the reproduction. Reset timer
                        t0 = t;
                        counter = 0;
                        disp('Underrun')
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
                        done(k) = isDone(obj.fileReader{k});
                    end
                    if all(done)
                        finish = true;
                        obj.playingState = playingStateClass('stopped');
                    end
                end
                
            end
            
        end
        
    end
    
    
    methods(Access = public)
        
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
        
        
        function obj = reproductor_plus()
            
            obj.playingState = playingStateClass('stopped');
            
            % Reading object
            obj.fileReader = {dsp.AudioFileReader()};
            
            % Processing object
            obj.processor = {processSignal('delayType', 'forward')};
            
            % Writing object
            obj.player = {audioPlayer};
            
            obj.comMatrix = 1;
            obj.setDefaultProperties();
            
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
                            release(obj); % It only executes if obj is locked. Make sure that always the state is playing, the obect is locked
                        case 'play'
                            % Two possibilities. Play the same song (reset)
                            % or play a new song
                            if ismember('fileName', p.UsingDefaults)
                                % Play current song from the beggining
                                % First, stop but without allowing the tuning
                                % of anything, i.e., don't release the objects
                                reset(obj);
                                % Then, it will play again since the state
                                % hasn't been changed
                            else
                                % Play a new audio file
                                % First, stop
                                release(obj);
                                % Then, assign the new file name
                                obj.audioFileName = p.Results.fileName;
                                % Then, play again
                                obj.playingState = playingStateClass('playing');
                                setup(obj, {}, {});
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
                            setup(obj, {}, {});
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
                            if ~isempty(obj.audioFileName)
                                setup(obj, {}, {});
                                obj.playingState = playingStateClass('playing');
                                obj.reproduce();
                            end
                        case 'resume'
                            % The same as play if there was no song
                            % specified
                            setup(obj, {}, {});
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
        
        function setProps(obj, parameters, values, indices )
            if obj.playingState == playingStateClass('stopped')
                
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
                            
                            % Set the value on the correspondent player and do
                            % the setup and release
                            player_this = obj.player{index};
                            player_this.driver = value;
                            setup(player_this, []);
                            release(player_this);
                            
                            % Get devices for the new driver and set it
                            audioDevices = player_this.getAvailableDevices();
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
                        case 'getDelayFun'                            
                            ind = obj.sub2LinkAbs(index(1), index(2));
                            obj.(parameter){ind} = value;
                        case 'getAttenFun'
                            ind = obj.sub2LinkAbs(index(1), index(2));
                            obj.(parameter){ind} = value;
                        case 'comMatrix'
                            obj.setCommutationMatrix(value);
                        case 'frameSizeReading'
                            obj.(parameter)(index) = value;
                        case 'Fs_player'
                            obj.(parameter)(index) = value;
                        case 'frameDuration'
                            obj.(parameter) = value;
                        otherwise
                            error('reproductor_plus:setProps', 'The first argument must be the name of an existing parameter')
                    end
                    
                end
            end
        end
        
        function reproduceNoRealTime(obj, t, delay, attenuation)
            if ~strcmp(obj.playingState, 'stopped'),  return;  end
            setup(obj, [], []);
            obj.playingState = playingStateClass('playing');
            
            Fs = obj.fileReader.SampleRate;
            numChann = obj.numChannels;
            numSamp = obj.frameSize;
            
            % Calculate samples where delay and attenuation changes
            t_Samp = ceil(t*Fs) + 1; % Samples were the position should change
            [t_Samp, ind] = unique(t_Samp); % Eliminate redundant information
            delay = delay(ind, :);
            attenuation = attenuation(ind, :);
            
            countSamples = 0;
            while ~isDone(obj.fileReader)
                
                % Find out which samples of change we need
                t_ind = find(t_Samp >= (countSamples + 1) & t_Samp <= (countSamples + numSamp));
                keySamples = [t_Samp(t_ind) - countSamples; numSamp + 1];
                if keySamples(1) > 1
                    keySamples = [1; keySamples];
                    prevInd = find(t_Samp < countSamples + 1, 1, 'last');
                    t_ind = [prevInd; t_ind];
                end
                
                % Give the delays and attenuations the convenient format for the processor
                % object
                delays = zeros(numSamp, numChann);
                attenuations = zeros(numSamp, numChann);
                for k = 1:numel(keySamples) - 1
                    currInd = keySamples(k):keySamples(k+1)-1;
                    delays(currInd, :) = repmat(delay(t_ind(k), :), numel(currInd), 1);
                    attenuations(currInd, :) = repmat(attenuation(t_ind(k), :), numel(currInd), 1);
                end
                
                step(obj, delays, attenuations);
                
                countSamples = countSamples + numSamp;
            end
            
            release(obj);
        end
          
    end
    
end

