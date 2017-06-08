function [ U_r_1, U_r_2 ] = kirchhoffIntegral( k, r_surface, dS, U, grad_U, r_measure )

numPointsSurface = size(r_surface, 1);
area_dS = modVec(dS);
n = scaleVector(dS, 1./area_dS);

rel_r = repmat(r_measure, numPointsSurface, 1) - r_surface;

% Method A
dist = modVec(rel_r);
G = exp(-1i*k*dist)./dist;
grad_G = scaleVector(-rel_r, exp(-1i*k*dist)./dist.^2.*(-1i*k - 1./dist));
U_r_1 = 1/(4*pi)*sum( area_dS.* (U.*dot(n, grad_G, 2) - G.*dot(n, grad_U, 2)) );

% % B1 - A1
% A1 = G.*dot(n, grad_U, 2);
% B1 = U.*dot(n, grad_G, 2);

% Method B
G = GreenFunction(rel_r, k);
grad_G = gradGreenFunction(-rel_r, k);
A = G .* dot(n, grad_U, 2);
B = U .* dot(n, grad_G, 2);
U_r_2 = 1/(4*pi)*sum( area_dS.* (A - B) );

% % A2 - B2
% A2 = A;
% B2 = B;

% U_surface_virtual = 1/(4*pi)*dS.*U_surface.*...
%     (- 1i*k*dot(scaleVector(r_p, 1./modVec(r_p)) - scaleVector(r_r0, 1./modVec(r_r0)), n, 2) ...
%     - dot(scaleVector(r_p, modVec(r_p).^(-2)), n, 2) ...
%     + dot(scaleVector(r_r0, modVec(r_r0).^(-2)), n, 2));


end

