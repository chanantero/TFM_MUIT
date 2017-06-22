Fs = 100;
t = (0:1/Fs:10)';
f = 0.1;
C = 1 + 1i;

x = zeros(numel(t)*4, 1);
x(numel(t)+1:2*numel(t)) = abs(C)*cos(2*pi*f*t + angle(C));
x(3*numel(t)+1:4*numel(t)) = abs(C)*cos(2*pi*f*t + angle(C));
tplot = (0:numel(x)-1)/Fs;
plot(tplot, x)

x2Comp = zeros(numel(t)*4, 1);
x2Comp(numel(t)+1:2*numel(t)) = C*exp(1i*2*pi*f*t);
x2Comp(3*numel(t)+1:4*numel(t)) = C*exp(1i*2*pi*f*t);
x2 = real(x2Comp);
tplot = (0:numel(x)-1)/Fs;
plot(tplot, x2, tplot, x)



[iq, xnew] = real2IQ(x, Fs, f);

ax = axes(figure);
plot3(tplot, real(xnew), imag(xnew), tplot, real(x2Comp), imag(x2Comp))
ax.XGrid = 'on';
ax.XLabel.String = 'Time';

A = exp(-1i*2*pi*f/Fs*(0:numel(x)-1)');
iq2 = x2Comp.*A;
plot3(tplot, real(A), imag(A))
plot3(tplot, real(iq2), imag(iq2))

XIQ = fft(x2Comp);
plot(abs(XIQ))
