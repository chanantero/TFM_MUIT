function value = gradGreenFunction(x, k, dim)
% dim. Dimension with size 3 that corresponds to the x, y, and z
% coordinates of the points in x

if nargin < 3
    dim = 2;
end

ndim = ndims(x);

dist = sqrt(sum(x.^2, dim));
greenValues = exp(-1i*k*dist)./dist;
scaleFactor = greenValues.*(-1i*k - 1./dist)./dist;

expVec = ones(1, ndim); expVec(dim) = 3;
value = repmat(scaleFactor, expVec) .* x;

end