%% Hilo infinito como fuente secundaria de monopolo
% Monopole source at (-Ds, 0 0)
% Hilo infinito along z direction on x=0 and y=0 for now.

limZ = [-1000 1000]; % Limits of integration on the infinite source

xm = 1:1:500; % Point where field will be calculated
yl = 0; % Position of linear source
k = 2;

Ds = 1:10;

ds = @(Ds, y, z) sqrt(Ds.^2 + y.^2 + z.^2); % Distance from primary source to secondary source
dr = @(xm, y, z) sqrt(xm.^2 + y.^2 + z.^2); % Distance from point of the secondary source and point of measure
ds_ = @(z) ds(Ds, yl, z);
dr_ = @(xm, z) dr(xm, yl, z);
% f = @(x, z) Ds./ds_(z).*(1i*k + 1./ds_(z)).*exp(-1i*k*(dr_(x, z) + ds_(z)))./(dr_(x, z).*ds_(z)); % Contribution from the source at point z to the field of point x
f = @(Ds, xm, z) Ds./ds(Ds, yl, z).*(1i*k + 1./ds(Ds, yl, z)).*exp(-1i*k*(dr(xm, yl, z) + ds(Ds, yl, z)))./(dr(xm, yl, z).*ds(Ds, yl, z)); % Contribution from the source at point z to the field of point x

I = zeros(numel(xm), numel(Ds));
for d = 1:numel(Ds)
    fprintf('%d/%d\n', d, numel(Ds))
    for p = 1:numel(xm);
        fAux = @(z) f(Ds(d), xm(p), z);
        I(p, d) = integral(fAux, limZ(1), limZ(2));
    end
end

ax = axes(figure, 'NextPlot', 'Add');
plot(ax, xm, real(I))
plot(ax, xm, imag(I))
plot(ax, xm, abs(I))

aux = I./repmat(I(:, 1), 1, size(I, 2));
plot(ax, xm, abs(aux))



ax = axes(figure, 'NextPlot', 'Add');
Iteo = -pi*1i*besselh(0, 2, k*xm);
plot(ax, x, real(Iteo), x, real(I))
plot(ax, x, imag(Iteo), x, imag(I))

ax = axes(figure, 'NextPlot', 'Add');
aux = I./repmat(Iteo.', 1, size(I, 2));
plot(ax, xm, abs(aux))

ax = axes;
plot(ax,