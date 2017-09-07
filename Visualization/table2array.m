function [parameters, varargout] = table2array(T, parameterNames, varargin)
% Assume that the table T has the correct format of a matrix that has
% been vectorized
    parameters = cellfun(@(paramName) unique(T.(paramName)), parameterNames, 'UniformOutput', false);
    parameterLengths = cellfun(@(param) numel(param), parameters);
    
    dataNames = varargin;
    varargout = table2dataArray(T, parameterNames, parameterLengths, dataNames);
end

function mats = table2dataArray(T, parameterNames, parameterLengths, dataNames)
    T = sortrows(T, flip(parameterNames));
%     mat = false(parameterLengths);
%     mat(:) = T.(dataName);
mats = cell(size(dataNames));
for k = 1:numel(dataNames)
    if ~isscalar(parameterLengths)
        mats{k} = reshape(T.(dataNames{k}), parameterLengths);
    else
        mats{k} = T.(dataNames{k});
    end
end

end
