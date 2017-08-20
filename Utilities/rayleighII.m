function [ U_r ] = rayleighII( k, r_surface, dS, U, r )

numPointsSurface = size(r_surface, 1);
area_dS = modVec(dS);
n = scaleVector(dS, 1./area_dS);

r_rs = repmat(r, numPointsSurface, 1) - r_surface;
r_rs_mod = modVec(r_rs);

% G = exp(-1i*k*r_rs_mod)./r_rs_mod;
grad_G = scaleVector(-r_rs, exp(-1i*k*r_rs_mod)./r_rs_mod.^2.*(-1i*k - 1./r_rs_mod));

U_r = 2/(4*pi)*sum( area_dS.*U.*dot(n, grad_G, 2) );

end