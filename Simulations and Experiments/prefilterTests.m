f = 0:1500;
c = 340;
g1 = 5.8284;
g2 = 0.5;
wc = 2*pi*500;
s = 1i*2*pi*f;
H_IIR = g2*(sqrt(g1)./wc*s + 1)./(1./(wc*sqrt(g1))*s + 1);
ax = axes(figure, 'NextPlot', 'Add');
plot(ax, f, sqrt(f/c));
plot(ax, f, abs(H_IIR))
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = '|H|';
legend(ax, 'Ideal', 'IIR')
yyaxis(ax, 'right')
plot(ax, f, rad2deg(angle(H_IIR)))
ax.YLabel.String = '\angle{H} (º)';

x=[0 1 2 3 4 5 6];
y=[0 .8415 .9093 .1411 -.7568 -.9589 -.2794];
polynom = LagrangePolynomial( x, y );

xvec = 0:0.1:6;
yvec = polyval(sum, xvec);
plot(xvec, yvec)
hold on
scatter(x,y)