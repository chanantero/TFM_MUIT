function distances = calcDistances(pointsA, pointsB)
% pointsA. Mx3 matrix.
% pointsB. Nx3 matrix.
% distances. MxN matrix.

M = size(pointsA, 1);
N = size(pointsB, 1);

pointAmat = repmat(permute(pointsA, [1 3 2]), [1 N]);
pointBmat = repmat(permute(pointsB, [3 1 2]), [M 1]);

distances = sqrt(sum((pointAmat - pointBmat).^2, 3));
end