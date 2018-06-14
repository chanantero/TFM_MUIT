function [ax, varargout] = drawArray(vectors, data, indepDims, varargin)

[repVectors, repData] =...
    filterArrayForRepresentation(vectors, data, indepDims, varargin{:});

p = inputParser;
p.KeepUnmatched = true;
addParameter(p, 'surfOrPlot', 'plot');
parse(p, varargin{:})
surfOrPlot = p.Results.surfOrPlot;

numIndepDims = length(indepDims);
ax = axes(figure);
if numIndepDims == 1
    plot(ax, repVectors{1}, repData);
    varargout = {repVectors{1}, repData};
elseif numIndepDims == 2
    if strcmp(surfOrPlot, 'surf')
        x = repVectors{1};
        y = repVectors{2};
        surf(ax, x, y, repData');
        varargout = {repVectors{1}, repVectors{2}, repData};
    elseif strcmp(surfOrPlot, 'plot')
        plot(ax, repVectors{1}, repData);
        leg = cellstr(num2str(repVectors{2})');
        legend(ax, leg)
        varargout = {repVectors{1}, repData};
    end    
else
    warning('drawArray:wrongNumberIndepDims', 'The number of representable independent dimensions must be 1 or 2')
end


end
