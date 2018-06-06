function IR = createPureDeltaImpulseResponse(delays, values, fs, numSamples)
% function IR = createPureDeltaImpulseResponse(dist, )
% Input arguments:
% - delays. Vector of delays.
% - values. Vector of values. The same size as delays.
% - fs. Scalar. Sampling frequencies.

N = numel(delays);
indDelta = floor(delays*fs) + 1;

IR = zeros(N, numSamples);
for n = 1:N
    if values(n) ~= 0 && indDelta(n) <= numSamples
        IR(n, indDelta(n)) = values(n);
    end
end

end