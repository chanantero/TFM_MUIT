%% Hankel and Bessel functions
% Calculate the field produced by a linear infinite source. This is, the
% source won't be punctual but one-dimensional, placed along the z axis (x
% = 0 and y = 0). As it is infinite, the generated field will only depend
% on the radius, but it will be independent on the angle or height in
% cilyndrical components.
% In order to calculate the field, two ways are used. The first one is by
% iusing the position z as variable of integration. The second will be the
% angle.
limZ = [-1000 1000]; % Limits of integration on the infinite source

x = 0.01:0.1:50; % Point where field will be calculated. Distance from the source
k = 1;

d = @(x, z) sqrt(z.^2 + x.^2); % Distance from point of the source and point of measure
f1 = @(x, z) exp(-1i*k*d(x, z))./d(x, z); % Contribution from the source at point z to the field of point x
f2 = @(x, alpha) exp(-1i*k*x./cos(alpha))./cos(alpha);

I1 = zeros(numel(x), 1);
I2 = zeros(numel(x), 1);
for l = 1:numel(x);
    fAux = @(z) f1(x(l), z);
    I1(l) = integral(fAux, limZ(1), limZ(2));
    
    fAux = @(alpha) f2(x(l), alpha);
    I2(l) = integral(fAux, -pi/2+0.1, pi/2-0.1);
end

ax = axes(figure, 'NextPlot', 'Add');
plot(ax, x, real(I1))%, h, imag(I2));
plot(ax, x, imag(I1))%, h, imag(I2));
plot(ax, x, real(I2))
plot(ax, x, imag(I2))

j0 = @(x) sinc(x);
j1 = @(x) sin(x)./(x.^2) - cos(x)./x;
j2 = @(x) (3./(x.^3) - 1./x).*sin(x) - 3./x.^2.*cos(x);
plot(x, j0(x), x, j1(x), x, j2(x))

J0 = besselj(0,x);
J1 = besselj(1,x);
J2 = besselj(2,x);
plot(x, J0, x, J1, x, J2)

Y0 = bessely(0,x);
Y1 = bessely(1,x);
Y2 = bessely(2,x);
plot(x, Y0, x, Y1, x, Y2)

H = besselh(1,x);
plot(x, imag(H))
plot(x, real(H))


[maxtab, mintab] = peakdet(imag(I1), 0.01, x);
[maxtab, mintab] = peakdet(real(I1), 0.01, x);


[maxtab1, mintab1] = peakdet(j0(x), 0.01, x);
[maxtab2, mintab2] = peakdet(j1(x), 0.01, x);
[maxtab3, mintab3] = peakdet(j2(x), 0.01, x);
[maxtab4, mintab4] = peakdet(imag(H), 0.01, x);
[maxtab5, mintab5] = peakdet(real(H), 0.01, x);
[maxtabJ0, mintabJ0] = peakdet(J0, 0.01, x);
[maxtabJ1, mintabJ1] = peakdet(J1, 0.01, x);
[maxtabJ2, mintabJ2] = peakdet(J2, 0.01, x);
[maxtabY0, mintabY0] = peakdet(Y0, 0.01, x);
[maxtabY1, mintabY1] = peakdet(Y1, 0.01, x);
[maxtabY2, mintabY2] = peakdet(Y2, 0.01, x);


maxtab(2:end,2)./maxtab(1:end-1,2)
mintab(2:end,2)./mintab(1:end-1,2)

maxtab1(2:end,2)./maxtab1(1:end-1,2)
maxtab2(2:end,2)./maxtab2(1:end-1,2)
maxtab3(2:end,2)./maxtab3(1:end-1,2)
maxtab4(2:end,2)./maxtab4(1:end-1,2)
maxtab5(2:end,2)./maxtab5(1:end-1,2)
maxtabJ0(2:end,2)./maxtabJ0(1:end-1,2)
maxtabJ1(2:end,2)./maxtabJ1(1:end-1,2)
maxtabY0(2:end,2)./maxtabY0(1:end-1,2)
maxtabY1(2:end,2)./maxtabY1(1:end-1,2)


mintab1(2:end,2)./mintab1(1:end-1,2)
mintab2(2:end,2)./mintab2(1:end-1,2)
mintab3(2:end,2)./mintab3(1:end-1,2)
mintab4(2:end,2)./mintab4(1:end-1,2)
mintab5(2:end,2)./mintab5(1:end-1,2)
mintabJ0(2:end,2)./mintabJ0(1:end-1,2)
mintabJ1(2:end,2)./mintabJ1(1:end-1,2)
mintabY0(2:end,2)./mintabY0(1:end-1,2)
mintabY1(2:end,2)./mintabY1(1:end-1,2)


% Máximo imag(I1) se corresponde con mínimo J0 y viceversa
plot(x, -J0*pi, x, imag(I1))

% Máximo real(I1) se corresponde con mínimo Y0 y viceversa
plot(x, -Y0*pi, x, real(I1))

Iteo = -pi * (bessely(0,k*x) + 1i*besselj(0, k*x));
Iteo = -1i*pi*conj(besselh(0, 1, k*x));
plot(x, real(Iteo), x, real(I1))
plot(x, imag(Iteo), x, imag(I1))

% Conclusiones
Iteo = -pi*1i*besselh(0, 2, k*x);
ax = axes(figure, 'NextPlot', 'Add');
plot(ax, x, real(Iteo), x, real(I2))
plot(ax, x, imag(Iteo), x, imag(I2))

% Phase stationary method
k = 1;
x = 0.1:0.1:100;
Iteo = -pi*1i*besselh(0, 2, k*x);
approx = exp(-1i*k*x)./sqrt(x)*sqrt(2*pi/(1i*k));
ax = axes(figure, 'NextPlot', 'Add');
plot(ax, x, real(Iteo), x, real(approx))
plot(ax, x, imag(Iteo), '--', x, imag(approx), '--')

rel = approx./Iteo;
ax = axes(figure, 'NextPlot', 'Add');
plot(ax, x, 10*log10(abs(rel).^2 - 1))
plot(ax, x, 10*log10(rad2deg(abs(angle(rel)))))
xlabel('kr')
legend('Amplitude error', 'Phase error')
