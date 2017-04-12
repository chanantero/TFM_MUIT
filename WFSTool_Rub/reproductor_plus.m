classdef reproductor_plus < matlab.System
    
    % The user can set this properties when the object is not locked, i.e,
    % when playingState is stopped
    properties(Nontunable, SetAccess = private)
        audioFileName % Cell string array with as many elements as fileReaders
        driver % Cell string array with as many elements as players
        device % Cell string array with as many elements as players
        getDelayFun % Cell array of functions with as many elements as players
        getAttenFun % Cell array of functions with as many elements as players
        comMatrix % Commutation matrix. numReaders x numPlayers.
        frameSizeReading % Array with as many elements as fileReaders 
        frameSizeWriting % Array with as many elements as players
    end
       
    properties(DiscreteState)
        count
    end
    
    % The user cannot set this properties directly, but are useful information that
    % can be viewed by other objects
    properties(SetAccess = private, SetObservable)%, AbortSet)
        playingState
        numChannels % Current number of output channels. It will be set each time a writing device changes
        numReaders
        numPlayers
        numLinks
    end

    % Private properties
    properties(SetAccess = private)
        fileReader % Cell Array
        processor % Cell Array as non-null elements of comMatrix
        player % Cell Array
    end
    
    methods(Access = protected)
        
        function setupImpl(obj)
            numRead = obj.numReaders;
            numPlay = obj.numPlayers;            
            indCom = find(obj.comMatrix ~= 0);
            [indReader, indPlayer] = ind2sub([numRead, numPlay], indCom);
            
            % Set propertties for the fileReaders
            Fs_reader = zeros(numRead, 1);
            for k = 1:numRead
                obj.fileReader{k}.Filename = obj.audioFileName(k);
                obj.fileReader{k}.SamplesPerFrame = obj.frameSizeReading(k);
                Fs_reader(k) = obj.fileReader{k}.SampleRate;
            end
            
            % Set properties for the processors
            for k = 1:obj.numLinks
                obj.processor{k}.Fs = Fs_reader(indReader(k));
                obj.processor{k}.numChannels = obj.numChannels(indPlayer(k));
            end
            
            % Set properties for the players
            Fs_players = zeros(numPlay, 1);
            for k = 1:numPlay
                Fs_players(k) = find(obj.comMatrix(:, k), 1, 'first');
            end
            for k = 1:numPlay
                obj.player{k}.frameSize = obj.frameSizeWriting(k);
                obj.player{k}.driver = obj.driver(k);
                obj.player{k}.device = obj.device(k);
                obj.player{k}.Fs = Fs_players(k);
            end
        end
        
        function numUnderrun = stepImpl(obj, delays, attenuations)
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
                Fs_reader = obj.processor{indReader(k)}.Fs;
                Fs_player = obj.player{indPlayer(k)}.Fs; 
                [P, Q] = rat(Fs_player/Fs_reader);
                audioProcFs = resample(audioProc{k}, P, Q);
            end
            
            % Sum outputs according to the commutation matrix
            % All processed signals corresponding to the same player must
            % have the same number of samples
            audioOutput = cell(numPlayer, 1);
            for k = 1:numPlayer
                ind = find(indPlayer == k);
                for l = 1:numel(ind)
                    audioOutput{k} = audioOutput{k} + obj.comMatrix(indCom(ind(l)))*audioProcFs{ind(l)};
                end
            end
            
            % Write to audio device buffer
            numUnderrun = zeros(numPlayer, 1);
            for k = 1:numPlayer
                numUnderrun(k) = step(obj.player{k}, audioOutput{k});
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
            % Check that the audioFileName exists
            a = dir(obj.audioFileName);
            assert(numel(a) > 0, 'reproductor:wrongProperty', 'The property audioFileName must specify an audio file that exists')
            
            % Check that the audioFileName has the right extension
            [~, ~, ext] = fileparts(obj.audioFileName);
            extensions = {'.wav', '.mp3'};
            assert(ismember(ext, extensions), 'reproductor:wrongProperty', ['The property audioFileName must have one of the next extensions: ', strjoin(extensions, ', ')])
            
        end
    end
    
    methods % Getters and setters
        function numPlayers = get.numPlayers(obj)
            numPlayers = numel(obj.players);
        end
        
        function numReaders = get.numReaders(obj)
            numReaders = numel(obj.fileReaders);
        end
        
        function numLinks = get.numLinks(obj)
            numLinks = sum(obj.comMatrix ~= 0);
        end
        
    end
    
    methods(Access = private)
        
        function setPropertiesFileReaders(obj, index)            
            if nargin == 1
                index = 1:obj.numReaders;
            end
            
            for k = index
                obj.fileReader{k}.Filename = obj.audioFileName(k);
                obj.fileReader{k}.SamplesPerFrame = obj.frameSizeReading(k);
            end
        end
        
        function setPropertiesProcessors(obj, index)
            if nargin == 1
                index = 1:obj.numLinks;
            end
            
            for k = 1:index
                obj.processor{k}.Fs = Fs_reader(indReader(k));
                obj.processor{k}.numChannels = obj.numChannels(indPlayer(k));
            end
        end
        
        function updateNumOutputChannels(obj, index)
            % The number of channels is the maximum number of
            % channels of the writing device
            if nargin == 1
                index = 1:obj.numPlayers;
            end
                        
            numChann = zeros(numel(index), 1);
            for k = index
                numChann(k) = obj.players{k}.numChannels;
            end
            
            obj.numChannels(index) = numChann;
        end
        
        function stop(obj)
            for k = 1:obj.numReaders
                release(obj.fileReaders{k})
            end
            
            for k = 1:obj.numLinks
                release(obj.processors{k})
            end
            
            for k = 1:obj.numPlayers
                release(obj.players{k})
            end
           
            obj.playingState = playingStateClass('stopped');
        end
        
        function reset_priv(obj)
            for k = 1:obj.numReaders
                reset(obj.fileReaders{k})
            end
            
            for k = 1:obj.numLinks
                reset(obj.processors{k})
            end
            
            for k = 1:obj.numPlayers
                reset(obj.players{k})
            end
            
            count = 0;
        end
        
        function setDefaultProperties(obj)
            for k = 1:obj.numReaders
                obj.audioFileName{k} = '';
                obj.frameSizeReading{k} = 1024*10;
            end
            
            for k = 1:obj.numWriters
                obj.frameSizeWriting{k} = 1024*10;
                obj.driver{k} = 'DirectSound';
                obj.device{k} = 'Default'; 
            end
            
            obj.updateNumOutputChannels();
        end
        
        function reproduce(obj)                
            % Timing control
            margin = 0.01;
            counter = 0;
            minPause = 0.01;
            minBufferDepletionTime = 0;
            t0 = tic;
            
            finish = false;
            while ~finish % Only reproduce if it is playing
                t = tic;
                
                if toc(t0) >= minBufferDepletionTime
                    delay = obj.getDelayFun();
                    attenuation = obj.getAttenFun();
                    delayMat = repmat(delay', obj.frameSizeReading, 1);
                    attenuationMat = repmat(attenuation', obj.frameSizeReading, 1);
                    
                    try
                        numUnderrun = step(obj, delayMat, attenuationMat);
                    catch
                        warning('There was some error with the step function of reproductor')
                        order.action = 'stop';
                        executeOrder(obj, order);
                        return;
                    end
                    
                    % Timing control
                    if numUnderrun > 0
                        % Interruption in the reproduction. Reset timer
                        t0 = t;
                        counter = 0;
                    end
                                           
                    counter = counter + 1;
                    minBufferDepletionTime = counter*obj.frameSizeReading/obj.player.Fs - margin;
                    pause(max(minPause, minBufferDepletionTime - toc(t0)));
%                     fprintf(max(minPause, minBufferDepletionTime - toc(t0)))
                else
                    pause(minPause);
                end
                
                % Finish condition
                if ~strcmp(obj.playingState, 'playing')
                    finish = true;
                elseif isDone(obj.fileReader)
                    finish = true;
                    obj.playingState = playingStateClass('stopped');
                end
                
            end
            
        end
        
    end
        
    methods(Access = public)
        
        function obj = reproductor_plus()
            
            obj.playingState = playingStateClass('stopped');
                        
            % Reading object
            obj.fileReaders = dsp.AudioFileReader();
            
            % Processing object
            obj.processors = processSignal('delayType', 'forward');
            
            % Writing object
            obj.players = audioPlayer;
            
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
                            release(obj); % It only executes if obj is locked. Make sure that always that the state is playing, the obect is locked
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
                                setup(obj, [], []);
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
                            setup(obj, [], []);
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
                                setup(obj, [], []);
                                obj.playingState = playingStateClass('playing');
                                obj.reproduce();
                            end
                        case 'resume'
                            % The same as play if there was no song
                            % specified
                            setup(obj, [], []);
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
        
        function setProps(obj, parameter, value, index )
            if obj.playingState == playingStateClass('stopped')
                
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
                        
                        % Get devices for the new driver
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
                        obj.(parameter){index} = value;
                    case 'getAttenFun'
                        obj.(parameter){index} = value;
                    case 'comMatrix'
                        obj.(parameter){index} = value;
                    case 'frameSizeReading'
                        obj.(parameter){index} = value;
                    case 'frameSizeWriting'
                        obj.(parameter){index} = value;
                    otherwise
                        error('reproductor_plus:setProps', 'The first argument must be the name of an existing parameter')
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
        
%         function setFrameSizeReading(obj, frameSize, index)
%             if nargin < 3
%                 index = 1;
%             end
%             
%             obj.setProps('frameSizeReading', frameSize, index);
%         end
%         
%         function setFrameSizeWriting(obj, frameSize, index)
%             if nargin < 3
%                 index = 1;
%             end
%             
%             obj.setProps('frameSizeWriting', frameSize, index);
%         end
%         
%         function set_getDelayFun(obj, getDelayFun, index)
%             if nargin < 3
%                 index = 1;
%             end
%             
%             obj.setProps('getDelayFun', getDelayFun, index);
%         end
%         
%         function set_getAttenFun(obj, getAttenFun)
%             if nargin < 3
%                 index = 1;
%             end
%             
%             obj.setProps('getAttenFun', getAttenFun, index);
%         end
%         
%         function devices = getWritingDevices(obj)
%             aux = audioDeviceWriter('Driver', obj.driver);
%             devices = getAudioDevices(aux);
%             % devices = obj.player.getAvailableDevices();
%         end
        
        
    end
    
end

