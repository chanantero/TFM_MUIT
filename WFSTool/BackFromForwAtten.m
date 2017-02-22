function backAtten = BackFromForwAtten(prevForwDelay, forwDelay, prevForwAtten, attenuation, storedBackAtten, Fs)
% We assume that the delay follows the theory in...

numPrevDelay = numel(prevForwDelay);
numBackDelay = numel(storedBackAtten);

firstSample = numBackDelay + 1;
lastSample = numel(forwDelay) + floor(forwDelay(end)*Fs); % With the floor(...) we make sure each interpolated point is calculated with one point at each side

backAtten = forward2BackwardDelay( [prevForwDelay; forwDelay], [prevForwAtten; attenuation], Fs, firstSample + numPrevDelay, lastSample + numPrevDelay );

backAtten = [storedBackAtten; backAtten'];

end
