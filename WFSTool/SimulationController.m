classdef SimulationController < handle
    % This class uses the WFSToolSimple class to make simulations easier to
    % perform.
    % Only one noise source at one frequency is accepted.
    % Only the WFS array of 96 loudspeakers is accepted.
    % The class can
    % perform various types of optimization in order to achieve
    % cancellation.
      
    properties
        cancelResults
        % cancelResults is a structure that describes a scenario, the
        % parameters that describe the cancellation attempts, and the
        % results of each cancellation attempt.
        % N. Number of cancellations attempts.
        % Basic variables that describe the scenario:
        % - NSRcoef
        % - NSVcoef
        % - NSRposition
        % - NSVposition
        % - WFSposition
        % - recPosition
        % - Frequency
        % Parameters of the cancellation:
        % - Type. Cell vector of N elements. Each element can have the next values:
        %   - 'Normal'. WFS coefficients created from the real noise source
        %   parameters and optimization.
        %   - 'OptVirtualNS'. The cancellation has been performed with the
        %   WFS coefficients generated from the virtual noise source
        %   parameters that minimizes the sum of the power of the signals
        %   received.
        % - OptimizationOptions. It only needs to exist when Type is
        % 'Normal'. Its a structure with 5 fields. Each field contains a
        % vector with as many elements as cancellation attempts
        % Revisar la estructura cancel Results, hay algo que me da mala
        % espina de la arquitectura que estoy eligiendo. Hay una razón para
        % decir esto, pero no sé articularla bien.
        % Results of the cancellation
        % - WFScoef
        % - recCoef
        % - recNScoef
        % - recWFScoef
    end
    
    properties(Dependent)
        % Microphones
        microPos % Microphone positions
        microCoef % It is the simulated field, but for interface purposes, it is included as a microphone variable
        microCoefNS % Simulated field with only the noise source (real, obviously) contributions
        microCoefWFS % Simulated field with only the WFS array contributions
        numMicro % Number of microphones
        
        % Noise Source
        amplitude
        amplitudeR
        amplitudeV
        phase
        phaseR 
        phaseV
        NScoef
        NSRcoef
        NSVcoef
        frequency % Frequency of the noise source
        channel % Channel of the audio device that will be used as a noise source
        NSposition
        NSRposition % Assumed real position of the noise source.
        NSVposition
        numNS
        
        % WFS
        WFSposition
        WFScoef
        numWFS
        
        % Cancellation attempts
        numCancellationAttempts
        
        % Visualization
        ax
        
        % Other
        Fs
    end
    
    properties(Access = public) % private
        WFSToolObj
    end
    
    % Getters and setters
    methods
        
        % Microphone variables
        function set.microPos(obj, value)
            numReceivers = size(value, 1);
            obj.WFSToolObj.setNumReceivers(numReceivers);
            obj.WFSToolObj.receiverPosition = value;
            obj.WFSToolObj.updateRecordPanelBasedOnVariables();
        end
        
        function microPos = get.microPos(obj)
            microPos = obj.WFSToolObj.receiverPosition;
        end
        
        function microCoef = get.microCoef(obj)
            microCoef = sum(obj.WFSToolObj.simulField, 2);
        end
        
        function microCoefNS = get.microCoefNS(obj)
            microCoefNS = obj.WFSToolObj.simulField(:, 1);
        end
        
        function microCoefWFS = get.microCoefWFS(obj)
            microCoefWFS = obj.WFSToolObj.simulField(:, 2);
        end
        
        function numMicro = get.numMicro(obj)
            numMicro = obj.WFSToolObj.numReceivers;
        end
        
        % Noise source variables
        function amplitude = get.amplitude(obj)
            amplitude = obj.WFSToolObj.amplitude;
        end
        
        function set.amplitude(obj, value)
            if size(value, 1) == 1
                obj.WFSToolObj.amplitude = [value; value];
            elseif size(value, 1) == 2
                obj.WFSToolObj.amplitude = value;
            end
            obj.WFSToolObj.updateReprodPanelBasedOnVariables();
        end
        
        function amplitudeR = get.amplitudeR(obj)
            amplitudeR = obj.WFSToolObj.amplitude(1);
        end
        
        function set.amplitudeR(obj, value)
            obj.WFSToolObj.amplitude(1) = value;
        end
        
        function amplitudeV = get.amplitudeV(obj)
            amplitudeV = obj.WFSToolObj.amplitude(2);
        end
        
        function set.amplitudeV(obj, value)
            obj.WFSToolObj.amplitude(2) = value;
        end
        
        function phase = get.phase(obj)
            phase = obj.WFSToolObj.phase;
        end
        
        function set.phase(obj, value)
            if size(value, 1) == 1
                obj.WFSToolObj.phase = [value; value];
            elseif size(value, 1) == 2
                obj.WFSToolObj.phase = value;
            end            
            obj.WFSToolObj.updateReprodPanelBasedOnVariables();
        end

        function phaseR = get.phaseR(obj)
            phaseR = obj.WFSToolObj.phase(1);
        end
        
        function set.phaseR(obj, value)
            obj.WFSToolObj.phase(1) = value;
        end
        
        function phaseV = get.phaseV(obj)
            phaseV = obj.WFSToolObj.phase(2);
        end
        
        function set.phaseV(obj, value)
            obj.WFSToolObj.phase(2) = value;
        end
        
        function NScoef = get.NScoef(obj)
            NScoef = obj.WFSToolObj.noiseSourceCoefficient;
        end
        
        function set.NScoef(obj, value)
            if size(value, 1) == 1
                obj.WFSToolObj.noiseSourceCoefficient = [value; value];
            else
                obj.WFSToolObj.noiseSourceCoefficient = value;
            end
        end
        
        function NSRcoef = get.NSRcoef(obj)
            NSRcoef = obj.WFSToolObj.noiseSourceCoefficient(1, :);
        end
        
        function set.NSRcoef(obj, value)
            obj.WFSToolObj.noiseSourceCoefficient(1, :) = value;
        end
        
        function NSVcoef = get.NSVcoef(obj)
            NSVcoef = obj.WFSToolObj.noiseSourceCoefficient(2, :);
        end
        
        function set.NSVcoef(obj, value)
            obj.WFSToolObj.noiseSourceCoefficient(2, :) = value;
        end
        
        function frequency = get.frequency(obj)
            frequency = obj.WFSToolObj.frequency(1);
        end
        
        function set.frequency(obj, value)
            obj.WFSToolObj.frequency = [value; value];
            obj.WFSToolObj.updateReprodPanelBasedOnVariables();
        end
        
        function NSposition = get.NSposition(obj)
            NSposition = obj.WFSToolObj.noiseSourcePosition;
        end
        
        function set.NSposition(obj, value)
            if size(value, 1) == 1
                obj.WFSToolObj.noiseSourcePosition = [value; value];
            elseif size(value, 1) == 2
                obj.WFSToolObj.noiseSourcePosition = value;
            end
        end
        
        function NSRposition = get.NSRposition(obj)
            NSRposition = obj.WFSToolObj.noiseSourcePosition(1, :);
        end
        
        function set.NSRposition(obj, value)
            obj.WFSToolObj.noiseSourcePosition(1, :) = value;
        end
        
        function NSVposition = get.NSVposition(obj)
            NSVposition = obj.WFSToolObj.noiseSourcePosition(2, :);
        end
        
        function set.NSVposition(obj, value)
            obj.WFSToolObj.noiseSourcePosition(2, :) = value;
        end
        
        function numNS = get.numNS(obj)
            numNS = obj.WFSToolObj.numNoiseSources;
        end
        
        % WFS
        function WFSposition = get.WFSposition(obj)
            WFSposition = obj.WFSToolObj.WFSarrayPosition;
        end
        
        function set.WFSposition(obj, value)
            numWFS = size(value, 1);
            if numWFS ~= obj.WFSToolObj.numSourcesWFSarray
                obj.WFSToolObj.setNumWFSarraySources(numWFS);
            end
            obj.WFSToolObj.WFSarrayPosition = value;
        end
        
        function WFScoef = get.WFScoef(obj)
            WFScoef = obj.WFSToolObj.WFSarrayCoefficient(:, 2);
        end
        
        function set.WFScoef(obj, value)
            obj.WFSToolObj.WFSarrayCoefficient(:, 2) = value;
        end
        
        function ax = get.ax(obj)
            ax = obj.WFSToolObj.ax;
        end
        
        function numCancellationAttempts = get.numCancellationAttempts(obj)
            numCancellationAttempts = numel(obj.cancelResults);
        end
        
        function numWFS = get.numWFS(obj)
            numWFS = obj.WFSToolObj.numSourcesWFSarray;
        end
        
        % Other
        function Fs = get.Fs(obj)
            Fs = obj.WFSToolObj.Fs;
        end
        
        function set.Fs(obj, value)
            obj.WFSToolObj.Fs = value;
        end
       
    end
    
    methods
        function obj = SimulationController()
            obj.setUp();
        end
        
        function setUp(obj)
            obj.WFSToolObj = WFSToolSimple;
            
            % Set noise source variables
            obj.WFSToolObj.setNumNoiseSources(2);
            obj.WFSToolObj.noiseSourceChannelMapping = [1; 0];
            obj.WFSToolObj.setVirtual([false; true]);
            obj.WFSToolObj.setReal([true; false]);
            
            % Set WFS array variables
            obj.WFSToolObj.setNumWFSarraySources(96);
            
            obj.WFSToolObj.updateReprodPanelBasedOnVariables();
            obj.WFSToolObj.updateRecordPanelBasedOnVariables();            
            
        end
        
        function setAcousticPaths(obj, varargin)
            if nargin == 2
               if strcmp(varargin{1}, 'experimental')
                   % Calculate experimental acoustic path
                   obj.WFSToolObj.reproduceAndRecordForAcousticPaths();
                   obj.WFSToolObj.calculateExperimentalAcousticPaths();
                   
                   obj.WFSToolObj.setLoudspeakerAcousticPath(obj.expAcPathStruct);
                   % The next step is necessary because the virtual source doesn't map onto
                   % any loudspeaker, so setLoudspeakerAcousticPaths leaves the acoustic paths
                   % of the virtual source with a value of 0.
                   obj.WFSToolObj.noiseSourceAcPathStruct.acousticPaths(:, 2) = obj.noiseSourceAcPathStruct.acousticPaths(:, 1); % The acoustic path of the virtual noise source is equal to the real one, for optimization purposes
               end
            else
                p = inputParser;

                addParameter(p, 'WFS', [])
                addParameter(p, 'NS', [])
                addParameter(p, 'loudspeakers', [])
                
                parse(p, varargin{:})
                
                for name_ = p.Parameters
                    name = name_{1};
                    if ~ismember(name, p.UsingDefaults)
                        switch name
                            case 'WFS'
                                x = p.Results.WFS;
                                if ischar(x) && strcmp(x, 'theoretical')
                                    obj.WFSToolObj.theoricWFSacousticPath();
                                else
                                    obj.WFSToolObj.WFSarrayAcPathStruct = x;
                                end                                
                            case 'NS'
                                x = p.Results.NS;
                                if ischar(x) && strcmp(x, 'theoretical')
                                    obj.WFSToolObj.theoricNoiseSourceAcousticPath();
                                else
                                    obj.WFSToolObj.noiseSourceAcPathStruct = x;
                                end              
                            case 'loudspeakers'
                                x = p.Results.loudspeakers;
                                obj.WFSToolObj.setLoudspeakerAcousticPath(x);
                        end
                    end
                end
            end
        end
        
        function cancel(obj, varargin)
            % The cancellation is performed in two steps.
            % 1) Get the theoric coefficients with the WFS
            % calculation.
            % 2) Optimize with those coefficients as a starting point.
            % Hence, one type of optimization is to just define the
            % optimization options of the 2nd step, and use the current
            % parameters to calculate the WFS coefficients in step 1.
            % The other type of optimization is a meta-optimization. This
            % is, we provide the optimization options of the 2nd step, and 
            % vary the parameters of the virtual noise source in the first
            % step to get different WFS coefficients.
            
            % A) Inputs
            noOptimisation = false;
            if nargin == 1
                numAttempts = 1;
                noOptimisation = true;
            else
                p = inputParser;
                
                addOptional(p, 'sourceFilter', [], @(x) all(ismember(x, {'NoFilter', 'Loudspeakers'})));
                addOptional(p, 'magConstraint', [], @(x) all(islogical(x)))
                addOptional(p, 'acousticPathType', [], @(x) all(ismember(x, {'Current', 'Theoretical'})));
                addOptional(p, 'grouping', [], @(x) all(ismember(x, {'Independent', 'AllTogether'})));
                addOptional(p, 'zerosFixed', false, @(x) all(islogical(x)));
                addParameter(p, 'testPoints', []);
                
                parse(p, varargin{:})
                
                sourceFilter = p.Results.sourceFilter;
                magConstraint = p.Results.magConstraint; % Maximum absolute value constraint (magnitude of loudspeaker complex coefficients less or equal than 1)
                acousticPathType = p.Results.acousticPathType;
                grouping = p.Results.grouping;
                zerosFixed = p.Results.zerosFixed;
                testPoints = p.Results.testPoints;
                defaultGrid = ismember('testPoints', p.UsingDefaults);
                
                count = [
                    numel(sourceFilter);
                    numel(magConstraint);
                    numel(acousticPathType);
                    numel(grouping);
                    numel(zerosFixed)];
                
                numAttempts = max(count);
                if any(count<numAttempts)
                    % Extend vectors
                    aux = {sourceFilter; magConstraint; acousticPathType; grouping; zerosFixed};
                    
                    for k = 1:5
                        if count(k) < numAttempts
                            % Extend vector
                            aux{k} = repmat(aux{k}(1), 1, 5);
                        end
                    end
                    
                    [sourceFilter, magConstraint, acousticPathType, grouping, zerosFixed] = aux{:};
                end
            end
            
            % B) Processing
            WFScoeff = zeros(obj.WFSToolObj.numSourcesWFSarray, obj.WFSToolObj.numNoiseSources, numAttempts);
            simulField = zeros(obj.WFSToolObj.numReceivers, obj.WFSToolObj.numNoiseSources, numAttempts);
            for k = 1:numAttempts
                % WFS cancellation
                if noOptimisation
                    obj.WFSToolObj.WFScalculation();
                else
                    if defaultGrid
                        obj.WFSToolObj.WFScalculation('SourceFilter', sourceFilter{k}, 'AcousticPath', acousticPathType{k}, 'Grouping', grouping{k}, 'maxAbsoluteValueConstraint', magConstraint(k), 'zerosFixed', zerosFixed(k));
                    else
                        obj.WFSToolObj.WFScalculation('SourceFilter', sourceFilter{k}, 'AcousticPath', acousticPathType{k}, 'Grouping', grouping{k}, 'maxAbsoluteValueConstraint', magConstraint(k), 'zerosFixed', zerosFixed(k), 'testPoints', testPoints);
                    end
                end
                
                % WFS array coefficients
                WFScoeff(:, :, k) = obj.WFSToolObj.WFSarrayCoefficient;
                                
                % Simulate
                obj.WFSToolObj.simulate();
                
                % Cancellation level (noise source 1 is real, noise source 2 is virtual)
                simulField(:, :, k) = obj.WFSToolObj.simulField;
            end
            
            WFScoeff = permute(WFScoeff(:, 2, :), [1, 3, 2]); % (numSourcesWFSarray x N). The WFS coefficients of the first frequency are all 0. It is reserved for the real noise source.
            WFSSimulField = permute(simulField(:, 2, :), [1, 3, 2]);
            totalSimulField = permute(sum(simulField, 2), [1 3 2]); % (numReceivers x N)
            
            % C) Format results
            S = repmat(obj.generateBasicExportStructure(), numAttempts, 1);
            for k = 1:numAttempts
                s = S(k);
                
                s.Type = 'normal';
                if noOptimisation
                    s.OptimizationOptions = 'No optimization';
                else
                    s.OptimizationOptions.sourceFilter = sourceFilter{k};
                    s.OptimizationOptions.magConstraint = magConstraint(k);
                    s.OptimizationOptions.acousticPathType = acousticPathType{k};
                    s.OptimizationOptions.grouping = grouping{k};
                    s.OptimizationOptions.zerosFixed = zerosFixed(k);
                end
                s.WFScoef = WFScoeff(:, k);
                s.recCoef = totalSimulField(:, k);
                s.recWFScoef = WFSSimulField(:, k);
                
                S(k) = s;
            end
                       
            obj.cancelResults = [obj.cancelResults; S];
        end
        
        function findBestVirtualSourceParameters(obj)
            % Find the theoric parameters for the virtual noise source that, applying
            % WFS cancellation with theoric acoustic path and unified optimization of loudspeakers, minimize the
            % magnitude of the resulting field using the experimental acoustic path
            
            % Optimize
            % The objective function accepts parameters as input arguments. It
            % returns the squared absolute value of the resulting field
            objectiveFunction = @(parameters) sum(abs(sum(noiseSourceParam2Field(parameters), 2)).^2);
            x0 = [obj.NSRposition, obj.amplitudeR, obj.phaseR]; % Initial value of parameters
            [xOpt, ~] = fminunc(objectiveFunction, x0); % Optimize
            
            % Set the parameters
            optPosition = xOpt(1:3);
            optAmplitude = xOpt(4);
            optPhase = xOpt(5);
            obj.WFSToolObj.amplitude(2) = optAmplitude;
            obj.WFSToolObj.phase(2) = optPhase;
            obj.WFSToolObj.noiseSourcePosition(2, :) = optPosition;
            obj.WFSToolObj.updateReprodPanelBasedOnVariables();
            
            % Set WFS array coefficients according to the optimized parameters. Respect
            % the generated theoric coefficients, don't ajust them independently or
            % something, just scale them to fit the theoric acoustic paths.
            obj.WFSToolObj.WFScalculation(...
                'SourceFilter', 'Loudspeakers',...
                'AcousticPath', 'Theoretical',...
                'Grouping', 'AllTogether',...
                'maxAbsoluteValueConstraint', false);
            
            % Simulate with the experimental acoustic path to see what is the resulting
            % field with these optimized theoric parameters of the virtual source
            obj.WFSToolObj.simulate();
                      
            s = obj.generateBasicExportStructure();
            s.Type = 'OptVirtualNS';
                        
            obj.cancelResults = [obj.cancelResults; s];
            
            function field = noiseSourceParam2Field(virtualNSparam)
                
                virtualNSparam = num2cell(virtualNSparam);
                [x, y, z, amp, ph] = virtualNSparam{:};
                
                % Set virtual noise source (2nd source) parameters
                pos = [x, y, z];
                obj.WFSToolObj.amplitude(2) = amp;
                obj.WFSToolObj.phase(2) = ph;
                obj.WFSToolObj.noiseSourcePosition(2, :) = pos;
                
                % WFS calculation
                obj.WFSToolObj.WFScalculation(...
                    'SourceFilter', 'Loudspeakers',...
                    'AcousticPath', 'Theoretical',...
                    'Grouping', 'AllTogether',...
                    'maxAbsoluteValueConstraint', true);
                
                % Simulate
                obj.WFSToolObj.simulate();
                
                % Output
                field = obj.WFSToolObj.simulField;
                
            end
            
        end
        
        function [NSposOpt, fVal] = findBestNoiseSourcePosition(obj)
            % The acoustic paths of the WFS array and the transmitted
            % signals are given, so the virtual noise source parameters are
            % irrelevant, since we don't perform any WFS calculation.
            % We optimize the real noise source position in order to
            % minimize the produced field.
            % It's like the noise source is cancellating the WFS field
            % instead of the opposite (which is the usual case)
            
            objectiveFunction = @(NSpos) NSpos2GlobalCancellation(NSpos);
            NSpos0 = obj.NSRposition; % Initial value of parameters
            [NSposOpt, fVal] = fminunc(objectiveFunction, NSpos0); % Optimize
            
            % Set the position of the solution
            
            function globalCancellation = NSpos2GlobalCancellation(NSpos)

                % Assign new position for the noise source
                obj.NSposition = NSpos;

                % Calculate theoric acoustic path for the noise source
                obj.setAcousticPaths('NS', 'theoretical');

                % Make acoustic paths effective
                obj.WFSToolObj.prepareOptimization();

                % Scale properly
                % Scale noise source coefficient. Given how cancellation is
                % done, I suspect it is not the appropiate way of
                % performing the optimization, since the noise source
                % coefficient appears in the denominator as well as the
                % numerator. It's more simple to scale the loudspeaker
                % coefficient. BUT! It's not the cancellation the thing in
                % which we are interested, isn't it? So, by now, I will
                % keep doing it this way.
                A = [obj.WFSToolObj.WFSarrayAcousticPath(:, :, 2), obj.WFSToolObj.noiseSourceAcousticPath(:, :, 1)];
                x0 = [obj.WFScoef; obj.NSRcoef; 0];
                indNSReal = obj.numWFS + 1;
                groups = {indNSReal};
                y = zeros(obj.numMicro, 1);
                xScaled = solveLinearSystem(A, y, groups, x0);
                obj.NScoef = xScaled(indNSReal);

                % Simulate
                obj.WFSToolObj.prepareSimulation();
                obj.WFSToolObj.simulate();
                
                globalCancellation = sum(abs(obj.microCoef).^2)/sum(abs(obj.microCoefNS).^2);
                
            end
        end
        
        function experimentalChecking(obj, varargin)
           
            for a = 1:obj.numCancellationAttempts
                obj.cancelResults(a).recNScoefExp = obj.reproduce(a, 'Noise'); % Received signal pulse coefficient matrix of the noise
                obj.cancelResults(a).recWFScoefExp = obj.reproduce(a, 'WFS');
                obj.cancelResults(a).recCoefExp = obj.reproduce(a, 'Total');
            end
                        
        end
        
        function recPulseCoefMat = reproduce(obj, index, source)
            
            obj.WFSToolObj.WFSarrayCoefficient(:, 2) = obj.cancelResults.WFScoef(:, index);
            
            switch source
                case 'Noise'
                    obj.WFSToolObj.setVirtual([false; false]); % No virtual sources
                    obj.WFSToolObj.setReal([true; false]); % Only real source
                    obj.WFSToolObj.updateReprodPanelBasedOnVariables();
                    
                case 'WFS'
                    obj.WFSToolObj.setVirtual([false; true]); % No virtual sources
                    obj.WFSToolObj.setReal([false; false]); % Only real source
                    obj.WFSToolObj.updateReprodPanelBasedOnVariables();

                    
                case 'Total'
                    obj.WFSToolObj.setVirtual([false; true]);
                    obj.WFSToolObj.setReal([true; false]);
                    obj.WFSToolObj.updateReprodPanelBasedOnVariables();
            end
            
            obj.WFSToolObj.WFScalculation(); % Update coefficents of loudspeakers. Only the channel of the real loudspeaker should have coefficient different than 0
            obj.WFSToolObj.reproduceAndRecord('main', 'soundTime', 2); % Simple reproduction of one pulse of 2 seconds
            
            sExp = obj.getExperimentalResultVariables();
            recPulseCoefMat = signal2pulseCoefficientMatrix(sExp.pulseLimits, obj.frequency, sum(sExp.pulseCoefMat, 3), sExp.recordedSignal, sExp.sampleRate).';
            
        end
        
        function s = generateBasicExportStructure(obj)
           s = SimulationController.generateExportStructure(...
                'NSRcoef', obj.NSRcoef,...
                'NSVcoef', obj.NSVcoef,...
                'WFScoef', obj.WFScoef,...
                'microCoef', obj.microCoef,...
                'microCoefNS', obj.microCoefNS,...
                'microCoefWFS', obj.microCoefWFS,...
                'NSRpos', obj.NSRposition,...
                'NSVpos', obj.NSVposition,...
                'WFSpos', obj.WFSposition,...
                'microPos', obj.microPos,...
                'Frequency', obj.frequency...
                );
        end
    end
    
    methods(Static)
        function s = generateExportStructure(varargin)
            p = inputParser;
            
            addParameter(p, 'NSRcoef', []);
            addParameter(p, 'NSVcoef', []);
            addParameter(p, 'WFScoef', []);
            addParameter(p, 'microCoef', []);
            addParameter(p, 'microCoefNS', []);
            addParameter(p, 'microCoefWFS', []);
            addParameter(p, 'NSRpos', []);
            addParameter(p, 'NSVpos', []);
            addParameter(p, 'WFSpos', []);
            addParameter(p, 'microPos', []);
            addParameter(p, 'Frequency', []);
            addParameter(p, 'OptimType', []);
            addParameter(p, 'OptimOptions', []);
            
            parse(p, varargin{:});
            
            s = struct;
            s.NSRcoef = p.Results.NSRcoef;
            s.NSVcoef = p.Results.NSVcoef;
            s.WFScoef = p.Results.WFScoef;
            s.recCoef = p.Results.microCoef;
            s.recNScoef = p.Results.microCoefNS;
            s.recWFScoef = p.Results.microCoefWFS;
            s.NSRposition = p.Results.NSRpos;
            s.NSVposition = p.Results.NSVpos;
            s.WFSposition = p.Results.WFSpos;
            s.recPosition = p.Results.microPos;
            s.Frequency = p.Results.Frequency;
            s.Type = p.Results.OptimType;
            s.OptimizationOptions = p.Results.OptimOptions;
        end
    end
    
end

