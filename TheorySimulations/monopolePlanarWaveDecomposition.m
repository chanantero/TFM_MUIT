% Direct sound enhancement by wave field synthesis, page 68.

k = 1;
z = 10;
G = @(x) exp(-1i*k*sqrt(z^2 + x.^2))./sqrt(z^2 + x.^2);

G_TF = @(x, kx) G(x).*exp(-1i*kx*x);

kx = -2:0.005:2;
TF = zeros(numel(kx), 1);
for l = 1:numel(kx)
    G_aux = @(x) G_TF(x, kx(l));
    TF(l) = integral(G_aux, -1000, 1000); 
end

% Theoretical value
kz = sqrt(k^2 - kx.^2); kz(imag(kz) ~= 0) = -kz(imag(kz) ~= 0);
TF_theo = sqrt(2*pi)./sqrt(kz * z).*exp(-1i*(kz*z + pi/4));

plot(kx, abs(TF), kx, abs(TF_theo))

plot(kx, angle(TF), kx, angle(TF_theo))

x = -1:0.01:1;
plot(x, 1./(1 - x.^4))

% Conclussion:
% The theoretical expression is valid for kz >> 1

