function [repVectors, repData] = filterArrayForRepresentation(vectors, data, indepDims, varargin)

p = inputParser;

p.KeepUnmatched = true;
addParameter(p, 'indepDimLimits', [])
addParameter(p, 'indepDimIndices', [])
addParameter(p, 'nonIndepDimValues', [])
addParameter(p, 'nonIndepDimIndices', [])

parse(p, varargin{:})

indepDimIndicesUser = ~ismember('indepDimIndices', p.UsingDefaults);
indepDimLimitsUser = ~ismember('indepDimLimits', p.UsingDefaults);
nonIndepDimIndicesUser = ~ismember('nonIndepDimIndices', p.UsingDefaults);
nonIndepDimValuesUser = ~ismember('nonIndepDimValues', p.UsingDefaults);

numDims = ndims(data);
numIndepDims = numel(indepDims);
numNonIndepDims = numDims - numIndepDims;
nonIndepDims = 1:numDims; nonIndepDims(ismember(nonIndepDims, indepDims)) = [];
sizeDims = size(data);

if indepDimIndicesUser
    indepDimIndices = p.Results.indepDimIndices;
else
    if indepDimLimitsUser
        indepDimLimits = p.Results.indepDimLimits;
        
        indepDimIndices = cell(numIndepDims, 1);
        for d = 1:numIndepDims
            vector = vectors{indepDims(d)};
            indepDimIndices{d} = vector >= indepDimLimits(d, 1) & vector <= indepDimLimits(d, 2);
        end
    else
        indepDimIndices = cell(numIndepDims, 1);
        for d = 1:numIndepDims
            indepDimIndices{d} = 1:sizeDims(indepDims(d));
        end
    end
end

if nonIndepDimIndicesUser
    nonIndepDimIndices = p.Results.nonIndepDimIndices;
else
    if nonIndepDimValuesUser
        nonIndepDimValues = p.Results.nonIndepDimValues;
        
        nonIndepDimIndices = zeros(numNonIndepDims, 1);
        for d = 1:numNonIndepDims
            vector = vectors{nonIndepDims(d)};
            objVal = nonIndepDimValues(d);
            if isfinite(objVal)
                [~, nonIndepDimIndices(d)] = min(abs(vector - objVal));
            else
                if objVal > 0
                    [~, nonIndepDimIndices(d)] = max(vector);
                else
                    [~, nonIndepDimIndices(d)] = min(vector);
                end
            end
        end
    else
        nonIndepDimIndices = ones(numNonIndepDims, 1);
    end
end

ind = cell(numDims, 1);
for d = 1:numIndepDims
    ind{indepDims(d)} = indepDimIndices{d};
end

for d = 1:numNonIndepDims
    ind{nonIndepDims(d)} = nonIndepDimIndices(d);
end

repData = data(ind{:}); % Filter the data
repData = permute(repData, [indepDims, nonIndepDims]); % Permute it
% The result is a matrix

repVectors = cell(numDims, 1);
for d = 1:numDims
    repVectors{d} = vectors{d}(ind{d});
end
repVectors = repVectors([indepDims, nonIndepDims]);

end