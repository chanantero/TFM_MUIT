function results = funcResult2grid( paramGrid, func, numOutputs )
% It takes a grid of points and applies a function R^N --> R^X, being X

% unknown
% numDim = numel(parameters);
% arraySize = cellfun(@(p) numel(p), parameters);
% 
% paramGrid = cell(numDim, 1);
% [paramGrid{:}] = ndgrid(parameters{:});

% We assume that the grid of every parameter has the same size. The i-th cell
% of paramGrid contains the values of the i-th parameter.
if nargin < 3
    numOutputs = 1;
end

numPoints = numel(paramGrid{1});
arraySize = size(paramGrid{1});

vecParamGrid = cellfun(@(p) num2cell(p(:)), paramGrid, 'UniformOutput', false);
vecParamGrid = cat(2, vecParamGrid{:});

if iscell(func)
    results = cell(numel(func), 1);
    for out = 1:numel(func)
        resultsCurr = cell(numOutputs, 1);
        results{out} = cell(numOutputs, 1);
        for k = 1:numOutputs
            results{out}{k} = cell(arraySize);
        end
        
        for k = 1:numPoints
            [resultsCurr{:}] = func{out}(vecParamGrid{k,:});
            for l = 1:numOutputs
                results{out}{l}(k) = resultsCurr(l);
            end
        end
    end
else
    resultsCurr = cell(numOutputs, 1);
    results = cell(numOutputs, 1);
    for k = 1:numOutputs
        results{k} = cell(arraySize);
    end
    
    for k = 1:numPoints
        [resultsCurr{:}] = func(vecParamGrid{k,:});
        for l = 1:numOutputs
            results{l}(k) = resultsCurr(l);
        end
    end
end
    

end

