pitch = 440;
SampleRate = 44100;
FrameSize = 44100;
A = 1 + 2*1i;
B = 65 - 0.26*1i;

t = (0:FrameSize-1)'/SampleRate;
x = real(A*exp(1i*2*pi*pitch*t));
x = [x; real(B*exp(1i*2*pi*pitch*t))];

windowSize = 30000;
values = complexCoefficient(x, SampleRate, pitch, windowSize);
t = (0:FrameSize*2 - 1)'/SampleRate;
plot(t(windowSize:end), real(values), t(windowSize:end), imag(values))
