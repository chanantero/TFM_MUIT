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
    end
    
    properties(Constant)
        numWFSsources = 96
    end
    
    properties(Dependent)
        % Microphones
        microPos % Microphone positions
        
        % Noise Source
        amplitude % Amplitude of the coefficient of the noise source
        phase % Phase of the coefficient of the noise source
        frequency % Frequency of the noise source
        channel % Channel of the audio device that will be used as a noise source
        NSposition % Assumed real position of the noise source.
        
        % Visualization
        ax
    end
    
    properties(Access = private)
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
        
        % Noise source variables
        function amplitude = get.amplitude(obj)
            amplitude = obj.WFSToolObj.amplitude(1);
        end
        
        function set.amplitude(obj, value)
            obj.WFSToolObj.amplitude = [value; value];
            obj.WFSToolObj.updateReprodPanelBasedOnVariables();
        end
        
        function phase = get.phase(obj)
            phase = obj.WFSToolObj.phase(1);
        end
        
        function set.phase(obj, value)
            obj.WFSToolObj.phase = [value; value];
            obj.WFSToolObj.updateReprodPanelBasedOnVariables();
        end

        function frequency = get.frequency(obj)
            frequency = obj.WFSToolObj.frequency(1);
        end
        
        function set.frequency(obj, value)
            obj.WFSToolObj.frequency = [value; value];
            obj.WFSToolObj.updateReprodPanelBasedOnVariables();
        end
        
        function NSposition = get.NSposition(obj)
            NSposition = obj.WFSToolObj.noiseSourcePosition(1, :);
        end
        
        function set.NSposition(obj, value)
            obj.WFSToolObj.noiseSourcePosition = [value; value];
        end
        
        function ax = get.ax(obj)
            ax = obj.WFSToolObj.ax;
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
            obj.WFSToolObj.setNumWFSarraySources(obj.numWFSsources);
            
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
            
            p = inputParser;
            
            addOptional(p, 'sourceFilter', [], @(x) all(ismember(x, {'NoFilter', 'Loudspeakers'})));
            addOptional(p, 'magConstraint', [], @(x) all(islogical(x)))
            addOptional(p, 'acousticPathType', [], @(x) all(ismember(x, {'Current', 'Theoretical'})));
            addOptional(p, 'grouping', [], @(x) all(ismember(x, {'Independent', 'AllTogether'})));
            addOptional(p, 'zerosFixed', [], @(x) all(islogical(x)));
            
            parse(p, varargin{:})
            
            sourceFilter = p.Results.sourceFilter;
            magConstraint = p.Results.magConstraint; % Maximum absolute value constraint (magnitude of loudspeaker complex coefficients less or equal than 1)
            acousticPathType = p.Results.acousticPathType;
            grouping = p.Results.grouping;
            zerosFixed = p.Results.zerosFixed;
            
            count = [
            numel(sourceFilter);
            numel(magConstraint);
            numel(acousticPathType);
            numel(grouping);
            numel(zerosFixed)];
        
            N = max(count);
            if any(count<N)
                % Extend vectors
                aux = {sourceFilter; magConstraint; acousticPathType; grouping; zerosFixed};
                
                for k = 1:5
                    if count(k) < N
                        % Extend vector
                        aux{k} = repmat(aux{k}(1), 1, 5);
                    end
                end
                
                [sourceFilter, magConstraint, acousticPathType, grouping, zerosFixed] = aux{:};
            end
            
            WFScoeff = zeros(obj.WFSToolObj.numSourcesWFSarray, obj.WFSToolObj.numNoiseSources, N);
            NScoeff = zeros(obj.WFSToolObj.numNoiseSources, obj.WFSToolObj.numNoiseSources, N);
            simulField = zeros(obj.WFSToolObj.numReceivers, obj.WFSToolObj.numNoiseSources, N);
            for k = 1:N
                % WFS cancellation
                obj.WFSToolObj.WFScalculation('SourceFilter', sourceFilter{k}, 'AcousticPath', acousticPathType{k}, 'Grouping', grouping{k}, 'maxAbsoluteValueConstraint', magConstraint(k), 'zerosFixed', zerosFixed(k));
                
                % WFS array coefficients
                WFScoeff(:, :, k) = obj.WFSToolObj.WFSarrayCoefficient;
                
                % NS coefficients
                NScoeff(:, :, k) = obj.WFSToolObj.noiseSourceCoefficient_complete;
                
                % Simulate
                obj.WFSToolObj.simulate();
                
                % Cancellation level (noise source 1 is real, noise source 2 is virtual)
                simulField(:, :, k) = obj.WFSToolObj.simulField;
            end
            
            WFScoeff = permute(WFScoeff(:, 2, :), [1, 3, 2]); % (numSourcesWFSarray x N). The WFS coefficients of the first frequency are all 0. It is reserved for the real noise source.
            NScoeff = diag(obj.WFSToolObj.noiseSourceCoefficient_complete);
            noiseSourceSimulField = simulField(:, 1, 1); % (numReceivers x 1). The real noise source coefficient doesn't change between optimizations, so we just select the first one
            WFSSimulField = permute(simulField(:, 2, :), [1, 3, 2]);
            totalSimulField = permute(sum(simulField, 2), [1 3 2]); % (numReceivers x N)
            
            s.WFScoef = WFScoeff;
            s.NScoef = NScoeff;
            s.recOnlyNoiseCoef = noiseSourceSimulField;
            s.recOnlyWFSCoef = WFSSimulField;
            s.recCoef = totalSimulField;
            
            obj.cancelResults = s;
        end
    end
    
end

