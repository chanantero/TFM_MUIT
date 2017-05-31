function U = sphericalWave_grad( k, CoefSources, rSources, r )

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
        difVector = repmat(permute(r, [1, 3, 2]), [1, numSources, 1])...
                  - repmat(permute(rSources, [3, 1, 2]), [numPoints, 1, 1]);
        
        distMatrix = sqrt(sum(difVector.^2, 3));
        
        CoefSourcesMat = repmat(permute(CoefSources, [2, 1]), [numPoints, 1]);
        A = CoefSourcesMat .* exp(-1i*k*distMatrix)./(distMatrix.^2).*(-1i*k - 1./distMatrix);
        
        U = permute(sum(repmat(A, [1 1 3]) .* difVector, 2), [1, 3, 2]);
        
    case 'loop'
        % With loop to avoid overflow of memory. It returns a column vector, not a
        % matrix
        U = zeros(numPoints, 3);
        for l = 1:numPoints
            % Calculate distance of every source to the point
            difVector = repmat(r(l,:), numSources, 1) - rSources;
            dist = sqrt(sum(difVector.^2, 2));
            
            A = CoefSources.*exp(-1i*k*dist)./(dist.^2).*(-1i*k - 1./dist);
            U(l, :) = sum(difVector * repmat(A, 1, 3), 1);
        end
end
end