function [ U, A, B ] = kirchhoffIntegral( k, r_surface, dS, Usurf, grad_Usurf, r_measure )

numMeasurePoints = size(r_measure, 1);
numPointsSurface = size(r_surface, 1);
area_dS = modVec(dS);
n = scaleVector(dS, 1./area_dS);

U = zeros(numMeasurePoints, 1);
for l = 1:numMeasurePoints
    r = r_measure(l, :);
    rel_r = repmat(r, numPointsSurface, 1) - r_surface;
    
    G = GreenFunction(rel_r, k);
    grad_G = gradGreenFunction(-rel_r, k); % The minus sign is because we must find the derivate of (x - x0)
    A = G .* sum(n.*grad_Usurf, 2);
    B = Usurf .* sum(n.*grad_G, 2);
    U(l) = 1/(4*pi)*sum( area_dS.* (A - B) );
end

end

