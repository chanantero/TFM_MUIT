function result = vecnorm(A, p, dim)
% It does the same as the official vecnorm function, but as it was
% implemented in Matlab2017b and I'm using Matlab2016a, I use this
% self-written version of the function

if nargin < 2
    p = 2;
end

if nargin < 3
    dim = 1;
end

result = (sum(A.^p, dim)).^(1/p);
end