G = @(k, r) exp(1i*k*sqrt(sum(r.^2, 2)))./sqrt(sum(r.^2, 2));

r0 = [2, 0, 0];

k = 2;
dx = 0.01;

derivNum = (G(k, r0+[dx, 0, 0]) - G(k, r0))/dx;


derivTeo = exp(1i*k*norm(r0))/norm(r0)^2*r0(1)*(1i*k - 1/norm(r0));


f = @(k,x) exp(1i*k*x)./x;

x = linspace(0,4,100);

y = f(k,x);

dxNum = diff(y)/(x(2)-x(1));
dxTeo = exp(1i*k*x)./x.*(1i*k - 1./x);

ax = axes(figure, 'NextPlot', 'Add');
plot(ax, x(1:end-1), real(dxNum))
plot(ax, x, real(dxTeo))


% Segunda derivada
k = 2;
r0 = [1e-1 0 0];
dx = 1e-6; dy = 1e-6; dz = 1e-6;

% Respecto a X
r = repmat(r0, 3, 1) + [[0; dx; 2*dx], zeros(3, 2)];
valuesX = G(k, r);
dx2 = (valuesX(3) - valuesX(2)*2 + valuesX(1))/dx^2;

% Respecto a Y
r = repmat(r0, 3, 1) + [zeros(3, 1), [0; dy; 2*dy], zeros(3, 1)];
valuesX = G(k, r);
dy2 = (valuesX(3) - valuesX(2)*2 + valuesX(1))/dy^2;

% Respecto a Z
r = repmat(r0, 3, 1) + [zeros(3, 2), [0; dz; 2*dz]];
valuesX = G(k, r);
dz2 = (valuesX(3) - valuesX(2)*2 + valuesX(1))/dz^2;

laplaciano = dx2 + dy2 + dz2;

laplaciano + k^2*G(k,r0) % == 0 en el caso ideal

dx2Teo = @(k, r) exp(1i*k*sqrt(sum(r.^2, 2))).*((1i*k - 1./sqrt(sum(r.^2, 2)))*(1./sum(r.^2, 2) + 1i*k*r(:,1).^2./sqrt(sum(r.^2, 2)).^3 - 2*r(:,1).^2./sqrt(sum(r.^2, 2)).^4) + r(:,1).^2./sqrt(sum(r.^2, 2)).^5);
dy2Teo = @(k, r) exp(1i*k*sqrt(sum(r.^2, 2))).*((1i*k - 1./sqrt(sum(r.^2, 2)))*(1./sum(r.^2, 2) + 1i*k*r(:,2).^2./sqrt(sum(r.^2, 2)).^3 - 2*r(:,2).^2./sqrt(sum(r.^2, 2)).^4) + r(:,2).^2./sqrt(sum(r.^2, 2)).^5);
dz2Teo = @(k, r) exp(1i*k*sqrt(sum(r.^2, 2))).*((1i*k - 1./sqrt(sum(r.^2, 2)))*(1./sum(r.^2, 2) + 1i*k*r(:,3).^2./sqrt(sum(r.^2, 2)).^3 - 2*r(:,3).^2./sqrt(sum(r.^2, 2)).^4) + r(:,3).^2./sqrt(sum(r.^2, 2)).^5);

laplacianoTeoCalc = dx2Teo(k,r0) + dy2Teo(k,r0) + dz2Teo(k,r0);
laplacianoTeoCalc + k^2*G(k,r0)

laplacianoTeo = @(k, r) -k^2*exp(1i*k*sqrt(sum(r.^2, 2)))./sqrt(sum(r.^2, 2));
laplacianoTeo(k,r0)

laplacianoTeo(k,r0) + k^2*G(k,r0)




secondDerivTeo = @(k, r) exp(1i*k*norm(r))./(norm(r).^2).*((1i*k-1./norm(r)).*(1i*k*r(:,1).^2./norm(r) - 2*r(:,1).^2./(norm(r).^2) + 1) + r(:,1).^2./norm(r).^3);
secondDerivNum = diff(y, 2)/(x(2)-x(1))^2;

secondDerivTeo(k, r0)


