function [originIndsAcum, destinationIndsAcum] = nestIndexing(varargin)
% The input arguments are pairs of indexes. Each pair consists of a vector
% of indexes that we will call originInds and another one called
% destinationInds. They represent a mapping between two vectors:
% B(destinationInds) = A(originInds);
% However, there is not only one pair, but N pairs. So, the relation is
% the next:
% A1(destinationInds{1}) = A0(originInds{1});
% A2(destinationInds{2}) = A1(originInds{2});
% ...
% Ai(destinationInds{i}) = Ai-1(originInds{i});
% ...
% AN(destinationInds{N}) = AN-1(originInds{N});
% We want to find the indexes that relate the last indices with the first
% ones, this is:
% AN(destinationIndsAcum) = A0(originIndsAcum);
% destinationIndsAcum and originIndsAcum won't have a number of elements
% bigger than any of the lengths of the individual index vectors.

if rem(nargin, 2) ~= 0
    error('The number of inputs must be even')
end

originInds = varargin(1:2:end);
destinationInds = varargin(2:2:end);
N = nargin/2; % Number of nested index relations
   
if N > 1
    % Recursion
    [flag, ind] = ismember(originInds{2}, destinationInds{1});
    destinationIndsAcum = destinationInds{2}(flag); % Acumulated destination indices
    originIndsAcum = originInds{1}(ind); % Acumulated origin indices
    aux = [originIndsAcum, originInds(3:end); destinationIndsAcum, destinationInds(3:end)];
    [originIndsAcum, destinationIndsAcum] = nestIndexing(aux{:});
else
     % Finish condition
    destinationIndsAcum = destinationInds{1};
    originIndsAcum = originInds{1};
end
    
end