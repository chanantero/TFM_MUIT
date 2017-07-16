function [iq, x_new] = real2IQ(x, SampleRate, pitchFreq, filterWidth)

[N, numChann] = size(x);
% FFT
X = fft(x);

clear('x');

% Delete the negative frequencies
X(ceil(N/2) + 1:end, :) = 0;
% Filter
if nargin == 4
    pitchFreqSamp = pitchFreq*N/SampleRate; % df = 1/duration = 1/(N/SampleRate)
    filterWidthSamp = filterWidth*N/SampleRate; % df = 1/duration = 1/(N/SampleRate)
    lowFreq = floor(pitchFreqSamp - filterWidthSamp);
    highFreq = ceil(pitchFreqSamp + filterWidthSamp);
    X(1:lowFreq-1, :) = 0;
    X(highFreq+1:end, :) = 0;
end
% Scale
X(2:end, :) = 2*X(2:end, :);
% Convert to time domain
x_new = ifft(X);

clear('X')

% Downconvert by the specified frequency
A = repmat(exp(-1i*2*pi*pitchFreq/SampleRate*(0:N-1)'), 1, numChann);
iq = x_new.*A;

end