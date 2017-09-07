function [ thresholds, changeDirections, allAboveThres, allUnderThres ] = searchThreshold( parameters, values, func, thresDim, threshold, tolerance )
% This function takes a grid of points in N dimensions, with an associated
% value to each point, with the function that gave as a result those
% values and, given a threshold value, finds via the shooting method the
% threshold along one of the dimensions for each of the points that form
% the rest N-1 dimensions.
% The inputs are:
% - parameters. Cell vector of N elements. Each one contains a vector of 
% points for each dimensions of the grid.
% - values. Array of N dimensions with the values of each point.
% - func. Function.
% - thresDim. Dimension on which we want to compute the threshold
% - threshold. Threshold value
% The outputs are:
% - thresholds. Array of N dimensions with the thresholds values. The
% thresDim-th dimensions has size 1.

numDim = numel(parameters);
sizeInputGrid = size(values);
sizeOutputGrid = sizeInputGrid; sizeOutputGrid(thresDim) = 1;
numPointsOutput = prod(sizeOutputGrid);
thresParam = parameters{thresDim};

thresholds = cell(sizeOutputGrid);
changeDirections = cell(sizeOutputGrid);
allAboveThres = false(sizeOutputGrid);
allUnderThres = false(sizeOutputGrid);
subs = cell(1, numDim);
for k = 1:numPointsOutput
    [subs{:}] = ind2sub(sizeOutputGrid, k);
    currParam = arrayfun(@(p) parameters{p}(subs{p}), 1:numDim, 'UniformOutput', false);
    currValues = values(subs{1:thresDim-1}, :, subs{thresDim+1:end});
    differ = currValues - threshold;
    signDiffer = diff(differ >= 0, 1);
    changeIntervals = find(signDiffer ~= 0);
    changeDirections{k} = signDiffer(changeIntervals);
    funcOneVariable = @(arg) func(currParam{1:thresDim-1}, arg, currParam{thresDim+1:numDim});
    thresholds{k} = zeros(1, numel(changeIntervals));
    for l = 1:numel(changeIntervals)     
        lim1 = thresParam(changeIntervals(l));
        lim2 = thresParam(changeIntervals(l)+1);
        limits = sort([lim1 lim2]);

        % Shooting method
        try
            thresholds{k}(l) = shootMethod(funcOneVariable, tolerance, limits, threshold);
        catch exception
            msg = ['Trying to execute shootMethod when k = ', num2str(k), ' and l = ', num2str(l), '.'];
            causeException = MException('searchThreshold:shootMethod', msg);
            exception = addCause(exception, causeException);
            rethrow(exception);
        end
    end
    
    if(isempty(changeIntervals))
        if(all(differ >= 0))
            allAboveThres(k) = true;
        elseif(all(differ < 0))
            allUnderThres(k) = true;
        else
            error('No se conoce la causa de que changeIntervals esté vacío')
        end    
    end
end

end

