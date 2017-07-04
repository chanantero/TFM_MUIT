function Mapprox = bestApproxMatrix( M, rank )
% Get the matrix with range == rank that has the minimum square error with the
% original matrix M (minimum square error = minimum frobenius norm)

[U, S, V] = svd(M);
Smod = zeros(size(S));
for k = 1:rank
    Smod(k, k) = S(k);
end
Mapprox = U*Smod*V;

end