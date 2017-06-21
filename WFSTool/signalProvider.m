classdef signalProvider < matlab.System
    
    % Public, tunable properties
    properties(Nontunable)
        mode % Read from audio file or generate sinusoidal signal
        SamplesPerFrame
        SampleRate % Sampling frequency. Not aplicable when mode == file
        
        % Variables that are exclusive for the mode == file
        FileName
        
        % Variables that are exclussive for the mode == sinusoidal
        amplitude
        phase
        frequency
        
        % Variables that are exclussive for the mode == custom
        customSignal
    end

    properties(SetAccess = private)
        stopSample % Sample at which the signal provider should reach the end. Only for mode == sinusoidal
    end
    
    properties(Dependent)
        numChannels
        numSamples
    end
    
    properties(DiscreteState)
        count
    end

    % Pre-computed constants
    properties(Access = private)
        fileReader
        providerFunction
    end

    methods(Access = protected)
        function setupImpl(obj)
            switch obj.mode
                case originType('file')
                    obj.providerFunction = @() obj.stepFileReader();
                    obj.fileReader.Filename = obj.FileName;
                    obj.fileReader.SamplesPerFrame = obj.SamplesPerFrame;
                    setup(obj.fileReader)
                case originType('sinusoidal')
                    obj.providerFunction = @() obj.toneGenerator();
                case originType('custom')
                    obj.providerFunction = @() obj.provideCustomSignal();
                otherwise
                    error('signalProvider:setup', 'The variable mode is not recognized')
            end
        end

        function y = stepImpl(obj)
            y = obj.providerFunction();
            obj.count = obj.count + 1;
        end

        function resetImpl(obj)
            obj.count = 0;
            reset(obj.fileReader);
        end
        
        function releaseImpl(obj)
            release(obj.fileReader);
        end
    end
    
    % Getters and setters
    methods
        function SampleRate = get.SampleRate(obj)
            if obj.mode == originType('file') || obj.mode == originType('custom') % Sinusoidal or custom
                SampleRate = obj.SampleRate;
            else %File
%                 SampleRate = obj.fileReader.SampleRate;
                try
                    inf = audioinfo(obj.FileName);
                    SampleRate = inf.SampleRate;
                catch
                    SampleRate = NaN;
                end
            end
        end
        
        function flag = isDone(obj)
            if obj.mode == originType('sinusoidal')
                if obj.count * obj.SamplesPerFrame >= obj.stopSample
                    flag = true;
                end
            elseif obj.mode == originType('custom')
                if obj.count * obj.SamplesPerFrame >= size(obj.customSignal, 1)
                    flag = true;
                else
                    flag = false;
                end
            elseif obj.mode == originType('file')
                flag = isDone(obj.fileReader);
            end
        end
        
        function numChannels = get.numChannels(obj)
            if obj.mode == originType('file')
                inf = audioinfo(obj.FileName);
                numChannels = inf.NumChannels;
            elseif obj.mode == originType('custom')
                numChannels = size(obj.customSignal, 2);
            else
                numChannels = 1;
            end
        end
        
        function numSamples = get.numSamples(obj)
            if obj.mode == originType('file')
                inf = audioinfo(obj.FileName);
                numSamples = inf.TotalSamples;
            elseif obj.mode == originType('custom')
                numSamples = size(obj.customSignal, 1);
            else
                numSamples = Inf;
            end
        end
        
%         function set.FileName(obj, fileName)
%             try
%                 obj.fileReader.Filename = fileName;
%             catch % e
% %                 warning('signalProvider:setFileName', e.message)         
%             end
%             obj.FileName = fileName;
%         end
    end
    
    methods
        function obj = signalProvider()
            obj.fileReader = dsp.AudioFileReader();
            obj.setDefaultProperties();
        end
    end
    
    methods(Access = private)
        function setDefaultProperties(obj)
            obj.mode = originType('file');
            obj.SamplesPerFrame = 44100;
            
            obj.FileName = '';
                        
            obj.amplitude = 1;
            obj.frequency = 1;
            obj.phase = 0;
            obj.SampleRate = 44100;      
        end
        
        function y = stepFileReader(obj)
            y = step(obj.fileReader);
        end
        
        function y = provideCustomSignal(obj)
            if obj.isDone()
                y = zeros(obj.SamplesPerFrame, obj.numChannels);
            else
                over = obj.count+1 * obj.SamplesPerFrame - obj.numSamples;
                if over > 0
                    ind = obj.count*obj.SamplesPerFrame + (1:obj.SamplesPerFrame - over)';
                    y = zeros(obj.SamplesPerFrame, obj.numChannels);
                    y(1:obj.SamplesPerFrame - over) = obj.customSignal(ind, :);                    
                else
                    ind = obj.count*obj.SamplesPerFrame + (1:obj.SamplesPerFrame)';
                    y = obj.customSignal(ind, :);
                end
            end
        end
        
        function y = toneGenerator(obj)
            if obj.isDone
                y = zeros(obj.SamplesPerFrame, 1);
            else
                A = obj.amplitude;
                Ph = obj.phase; % Initial phase
                f = obj.frequency;
                t = (obj.count*obj.SamplesPerFrame + (0:obj.SamplesPerFrame-1)')/obj.SampleRate;
                
                y = A*cos(2*pi*f*t + Ph);
            end
        end
        
    end
end
