classdef signalProvider < matlab.System
    
    % Public, tunable properties
    properties(Nontunable)
        mode % Read from audio file or generate sinusoidal signal
        SamplesPerFrame
        
        % Variables that are exclusive for the mode == file
        FileName
        
        % Variables that are exclussive for the mode == sinusoidal
        amplitude
        phase
        frequency
        SampleRate % Sampling frequency
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
                    obj.providerFunction = @() obj.signalGenerator();
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
            if obj.mode % Sinusoidal
                SampleRate = obj.SampleRate;
            else %File
                SampleRate = obj.fileReader.SampleRate;
            end
        end
        
        function flag = isDone(obj)
            if obj.mode
                flag = false;
            else
                flag = isDone(obj.fileReader);
            end
        end
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
        
        function y = signalGenerator(obj)
            A = obj.amplitude;
            Ph = obj.phase; % Initial phase
            f = obj.frequency;
            t = (obj.count*obj.SamplesPerFrame + (0:obj.SamplesPerFrame-1)')/obj.SampleRate;
            
            y = A*cos(2*pi*f*t + Ph);
        end
    end
end
