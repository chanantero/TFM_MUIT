function [ varargout ] = ndgrid_multClass( varargin )
% Make ndgrid but without worrying about the class of the vectors

indicesVec = cellfun(@(vec) 1:numel(vec), varargin, 'UniformOutput', false);

indicesMats = cell(numel(varargin), 1);
[indicesMats{:}] = ndgrid(indicesVec{:});

outputMats = cell(numel(varargin), 1);
for k = 1:numel(varargin)
    outputMats{k} = varargin{k}(indicesMats{k});
end

varargout = outputMats;

end

