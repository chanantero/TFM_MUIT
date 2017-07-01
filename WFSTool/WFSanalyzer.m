classdef WFSanalyzer
    
    properties
        ax
        fig
    end
    
    methods
        function obj = WFSanalyzer()
            obj.fig = figure;
        end
        
        function representExperimentalAndSimulated(obj, rat)
            
        end
        
        function analyzeSignal(obj, x, xSampleRate)
            clf(obj.fig);
            
            posAxSignal = [0 0.5 0.5 0.5];
            posAxSpectrum= [0.5 0.5 0.5 0.5];
            posAxIQ= [0 0 0.5 0.5];
            posAxDiff= [0.5 0 0.5 0.5];
            axSignal = axes(obj.fig, 'Units', 'normalized', 'OuterPosition', posAxSignal);
            axSpectrum = axes(obj.fig, 'Units', 'normalized', 'OuterPosition', posAxSpectrum);
            axIQ = axes(obj.fig, 'Units', 'normalized', 'OuterPosition', posAxIQ);
            axDiff = axes(obj.fig, 'Units', 'normalized', 'OuterPosition', posAxDiff);
            
            % Represent signal
            numSamplesX = size(x, 1);
            t = (0:numSamplesX-1)'/xSampleRate;            
            plot(axSignal, t, x);
           
            % Represent spectrum
            X = fft(x);
            df = (numSamplesX/xSampleRate)^(-1);
            f = (0:numSamplesX-1)*df;
            plot(axSpectrum, f, abs(X));
            
            % Represent IQ signal
            f = 600;
            iq = real2IQ(x, xSampleRate, f, 100);
            plot(axIQ, t, abs(iq))
            
            % Represent difference with recovered signal from IQ signal
            recovered = real(iq.*exp(1i*2*pi*f*t));
            difference = x - recovered;
            plot(axDiff, t, abs(difference))
        end
        
        function representRecordedSignal(obj, y, ySampleRate)            
%             axReprod = axes(obj.fig, 'Units', 'normalized', 'Position', [0.1 0.6 0.8 0.3]);
            axRecord = axes(obj.fig, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.3]);
            
%             numReprodSamples = size(x, 1);
%             tReprod = (0:numReprodSamples - 1)'/xSampleRate;
            
            numRecordSamples = size(y, 1);
            tRecord = (0:numRecordSamples - 1)'/ySampleRate;
            
%             plot(axReprod, tReprod, x)
            plot(axRecord, tRecord, y)
            
%             obj.ax = [axReprod; axRecord];
            obj.ax = axRecord;
        end
    end
    
end

