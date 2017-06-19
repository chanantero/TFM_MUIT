function iq = real2IQ(x, SampleRate, pitchFreq, windowSize)

if nargin == 3
    mode = 'accurate';
elseif nargin == 4
    mode = 'fast';
else
    error('real2IQ:wrongNumberOfInputs', 'Wrong number of inputs')
end

switch mode
    case 'fast'
        [N, numChann] = size(x);
        
        A = x .* repmat(cos(2*pi*pitchFreq*(0:N-1)'/SampleRate), 1, numChann);
        B = x .* repmat(sin(2*pi*pitchFreq*(0:N-1)'/SampleRate), 1, numChann);
        
        Acum = cumsum(A);
        Bcum = cumsum(B);
        realPart = (Acum(windowSize:end) - Acum(1:end-windowSize+1))*2/windowSize;
        imagPart = (Bcum(windowSize:end) - Bcum(1:end-windowSize+1))*(-2)/windowSize;
        
        iq = realPart + 1i*imagPart;
    case 'accurate'
        [N, numChann] = size(x);
        % FFT
        X = fft(x);
        % Delete the negative frequencies
        X(ceil(N/2) + 1:end, :) = 0;
        % Scale
        X(2:end) = 2*X(2:end);
        % Convert to time domain
        x_new = ifft(X);
        % Downconvert by the specified frequency
        A = repmat(exp(-1i*2*pi*pitchFreq/SampleRate*(0:N-1)'), 1, numChann);
        iq = x_new.*A;
end

end