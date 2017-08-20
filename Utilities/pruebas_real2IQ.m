pitch = 440;
SampleRate = 44100;
FrameSize = 44100;
A = 1 + 2*1i;
B = 6.5 - 0.26*1i;

t = (0:FrameSize-1)'/SampleRate;
x1 = zeros(numel(t), 1);
x2 = real(A*exp(1i*2*pi*pitch*t));
x3 = real(B*exp(1i*2*pi*pitch*t));
x = [x1; x2; x1; x3; x2];

windowSize = 30000;
values = real2IQ(x, SampleRate, pitch, windowSize);
t = (0:numel(x) - 1)'/SampleRate;
plot(t(windowSize:end), real(values), t(windowSize:end), imag(values))

values = real2IQ(x, SampleRate, pitch);
t = (0:numel(x) - 1)'/SampleRate;
plot(t, real(values), t, imag(values))
plot(real(values.*exp(1i*2*pi*pitch*t)))



