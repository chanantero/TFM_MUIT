function U = sphericalWave( k, CoefSources, rSources, r )

assert(size(rSources, 2) == 3, 'sphericalWave:wrongInput', 'size(rSources, 2) == 3')
assert(size(r, 2) == 3, 'sphericalWave:wrongInput', 'size(r, 2) == 3')

numPoints = size(r, 1);
numSources = size(rSources, 1);

if numPoints*numSources > 1e7
    mode = 'loop';
else
    mode = 'array';
end

switch mode
    case 'array'
        % With arrays
        distMatrix = zeros(numPoints, numSources);
        for l = 1:numSources
            difVector = r - repmat(rSources(l, :), numPoints, 1);
            distMatrix(:, l) = sqrt(sum(difVector.^2, 2));
        end
        
        U = repmat(CoefSources.', numPoints, 1).*exp(-1i*k*distMatrix)./distMatrix;
        U = sum(U, 2);
        
    case 'loop'
        % With loop to avoid overflow of memory. It returns a column vector, not a
        % matrix
        U = zeros(numPoints, 1);
        for l = 1:numPoints
            % Calculate distance of every source to the point
            difVector = repmat(r(l,:), numSources, 1) - rSources;
            dist = sqrt(sum(difVector.^2, 2));
            
            U(l) = CoefSources.'*(exp(-1i*k*dist)./dist);
        end
end
end