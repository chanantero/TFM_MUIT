function y = shiftAndCrop(x, shiftPos, newLength)

N = size(x, 1);

y = zeros(newLength, 1);

firstYind = max(1, 1 + shiftPos);
lastYind = min(newLength, N + shiftPos);

firstXind = max(1, 1 - shiftPos);
lastXind = min(N, newLength - shiftPos);

y(firstYind:lastYind) = x(firstXind:lastXind);

end
