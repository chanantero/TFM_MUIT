function vectors = scaleVector(vectors, scalingFactors)
% Multiplies N 3-dimensional vectors by N scalars

% Input arguments:
% - vectors. Nx3 array. Each row is a 3-element vector
% - scalingFactors. Nx1 array. The i-th element is the factor by which the
% i-th vector or 'vectors' is going to be multiplied
%
% Output arguments:
% - vectors. Scaled vectors

vectors = vectors.*repmat(scalingFactors, 1, 3);

end