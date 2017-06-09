function value = complexCoefficient(x, SampleRate, pitchFreq, windowSize)

[N, numChann] = size(x);

A = x .* repmat(cos(2*pi*pitchFreq*(0:N-1)'/SampleRate), 1, numChann);
B = x .* repmat(sin(2*pi*pitchFreq*(0:N-1)'/SampleRate), 1, numChann);

Acum = cumsum(A);
Bcum = cumsum(B);
realPart = (Acum(windowSize:end) - Acum(1:end-windowSize+1))*2/windowSize;
imagPart = (Bcum(windowSize:end) - Bcum(1:end-windowSize+1))*(-2)/windowSize;

value = realPart + 1i*imagPart;

end