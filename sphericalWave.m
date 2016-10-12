function U = sphericalWave( k, CoefSources, rSources, r )

assert(size(rSources, 2) == 3, 'sphericalWave:wrongInput', 'size(rSources, 2) == 3')
assert(size(r, 2) == 3, 'sphericalWave:wrongInput', 'size(r, 2) == 3')

numPoints = size(r, 1);
numSources = size(rSources, 1);

distMatrix = zeros(numPoints, numSources);
for k = 1:numSources
    difVector = r - repmat(rSources(k, :), numPoints, 1);
    distMatrix(:, k) = sqrt(sum(difVector.^2, 2));
end

U = repmat(CoefSources', numPoints, 1).*exp(-1i*k*distMatrix)./distMatrix;

end

