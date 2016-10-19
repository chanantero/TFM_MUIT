function [ U_r ] = kirchhoffIntegral( k, r_surface, dS, U, grad_U, r )

numPointsSurface = size(r_surface, 1);
area_dS = modVec(dS);
n = scaleVector(dS, 1./area_dS);

r_rs = repmat(r, numPointsSurface, 1) - r_surface;
r_rs_mod = modVec(r_rs);

G = exp(-1i*k*r_rs_mod)./r_rs_mod;
grad_G = scaleVector(-r_rs, exp(-1i*k*r_rs_mod)./r_rs_mod.^2.*(-1i*k - 1./r_rs_mod));

U_r = 1/(4*pi)*sum( area_dS.* (U.*dot(n, grad_G, 2) - G.*dot(n, grad_U, 2)) );

% U_surface_virtual = 1/(4*pi)*dS.*U_surface.*...
%     (- 1i*k*dot(scaleVector(r_p, 1./modVec(r_p)) - scaleVector(r_r0, 1./modVec(r_r0)), n, 2) ...
%     - dot(scaleVector(r_p, modVec(r_p).^(-2)), n, 2) ...
%     + dot(scaleVector(r_r0, modVec(r_r0).^(-2)), n, 2));


end

