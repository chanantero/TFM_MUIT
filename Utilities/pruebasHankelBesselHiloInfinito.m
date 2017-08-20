limX = [-1000 1000];

h = 0.01:0.1:50;
k = 1;

d = @(x, h) sqrt(h^2 + x.^2);

I1 = zeros(numel(h), 1);
% I2 = zeros(numel(h), 1);
for l = 1:numel(h);
    f1 = @(x) exp(-1i*k*d(x, h(l)))./d(x, h(l));
    I1(l) = integral(f1, limX(1), limX(2));
    
%     f2 = @(alpha) exp(-1i*k*h(l)./cos(alpha))./cos(alpha);
%     I2(l) = integral(f2, -pi/2+0.1, pi/2-0.1);
end

plot(h, real(I1))%, h, imag(I2));
plot(h, imag(I1))%, h, imag(I2));

j0 = @(x) sinc(x);
j1 = @(x) sin(x)./(x.^2) - cos(x)./x;
j2 = @(x) (3./(x.^3) - 1./x).*sin(x) - 3./x.^2.*cos(x);
plot(h, j0(h), h, j1(h), h, j2(h))

J0 = besselj(0,h);
J1 = besselj(1,h);
J2 = besselj(2,h);
plot(h, J0, h, J1, h, J2)

Y0 = bessely(0,h);
Y1 = bessely(1,h);
Y2 = bessely(2,h);
plot(h, Y0, h, Y1, h, Y2)

H = besselh(1,h);
plot(h, imag(H))
plot(h, real(H))


[maxtab, mintab] = peakdet(imag(I1), 0.01, h);
[maxtab, mintab] = peakdet(real(I1), 0.01, h);


[maxtab1, mintab1] = peakdet(j0(h), 0.01, h);
[maxtab2, mintab2] = peakdet(j1(h), 0.01, h);
[maxtab3, mintab3] = peakdet(j2(h), 0.01, h);
[maxtab4, mintab4] = peakdet(imag(H), 0.01, h);
[maxtab5, mintab5] = peakdet(real(H), 0.01, h);
[maxtabJ0, mintabJ0] = peakdet(J0, 0.01, h);
[maxtabJ1, mintabJ1] = peakdet(J1, 0.01, h);
[maxtabJ2, mintabJ2] = peakdet(J2, 0.01, h);
[maxtabY0, mintabY0] = peakdet(Y0, 0.01, h);
[maxtabY1, mintabY1] = peakdet(Y1, 0.01, h);
[maxtabY2, mintabY2] = peakdet(Y2, 0.01, h);


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
plot(h, -J0*pi, h, imag(I1))

% Máximo real(I1) se corresponde con mínimo Y0 y viceversa
plot(h, -Y0*pi, h, real(I1))

Iteo = -pi * (bessely(0,h) + 1i*besselj(0, h));
Iteo = -1i*pi*conj(besselh(0, 1, h));
plot(h, real(Iteo), h, real(I1))
plot(h, imag(Iteo), h, imag(I1))